use crate::ancestors::{AncestorArray, AncestralSequence, ANCESTRAL_STATE, DERIVED_STATE};
use crate::variants::{SequencePosition, VariantIndex, VariantSite};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::io::Write;
use std::path::Path;
use std::{io, mem};
use twox_hash::XxHash64;

/// Generates ancestral sequences for a given set of variant sites. The ancestral sequences are
/// generated using heuristic methods that use a small number of variant sites as focal sites and
/// infer the ancestral state for surrounding sites. For each set of focal sites, a single ancestral
/// sequence is generated.
pub struct AncestorGenerator {
    sites: Vec<VariantSite>,
}

impl AncestorGenerator {
    /// Create a new ancestor generator from an iterator over variant sites.
    pub fn from_iter(iter: impl Iterator<Item = VariantSite>) -> Self {
        Self {
            sites: iter.filter(|site| Self::is_valid_site(site)).collect(),
        }
    }

    /// Create a new ancestor generator without variant sites.
    pub fn empty() -> Self {
        Self { sites: Vec::new() }
    }

    /// Create a new ancestor generator without variant site, but with a given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            sites: Vec::with_capacity(capacity),
        }
    }

    /// Whether a site is valid for the generator. A site is valid if it is not a singleton and
    /// is biallelic.
    fn is_valid_site(site: &VariantSite) -> bool {
        !site.is_singleton && site.is_biallelic
    }

    /// Add a variant site to the generator.
    pub fn add_site(&mut self, site: VariantSite) {
        if Self::is_valid_site(&site) {
            self.sites.push(site);
        }
    }

    /// Extract the sequence positions of the variant sites.
    pub fn variant_positions(&self) -> Vec<SequencePosition> {
        self.sites.iter().map(|s| s.position).collect()
    }

    /// For a given set of focal sites, compute an ancestor that uses those focal sites. The focal
    /// site is given by a sorted set of indices into the set of variant sites.
    fn generate_ancestor(&self, focal_sites: &[usize]) -> AncestralSequence {
        debug_assert!(!focal_sites.is_empty());
        debug_assert!(focal_sites.windows(2).all(|sites| sites[0] < sites[1]));

        let mut ancestral_sequence = AncestralSequence::from_ancestral_state(
            self.sites.len(),
            self.sites[focal_sites[0]].relative_age,
        );
        let sites = self.sites.len();

        // extend ancestor to the left of the first focal site
        let modified_left = self.extend_ancestor(
            &mut self
                .sites
                .iter()
                .enumerate()
                .rev()
                .skip(sites - focal_sites[0]),
            focal_sites[0],
            &mut ancestral_sequence,
            true,
        );

        // infer ancestor between focal sites.
        for foc_index in 0..focal_sites.len() - 1 {
            let focal_site_i = focal_sites[foc_index];
            let focal_site_j = focal_sites[foc_index + 1];
            self.extend_ancestor(
                &mut self
                    .sites
                    .iter()
                    .enumerate()
                    .skip(focal_site_i + 1)
                    .take(focal_site_j - focal_site_i - 1),
                focal_site_i,
                &mut ancestral_sequence,
                false,
            );
        }

        // extend ancestor to the right of the last focal site
        let last_focal_site = *focal_sites.last().unwrap();
        let modified_right = self.extend_ancestor(
            &mut self.sites.iter().enumerate().skip(last_focal_site + 1),
            last_focal_site,
            &mut ancestral_sequence,
            true,
        );

        // TODO truncate ancestral sequence to save memory. this should be implemented using functionality from BitVec
        ancestral_sequence.focal_sites = focal_sites.iter().map(|&f| VariantIndex(f)).collect();
        ancestral_sequence.start = VariantIndex(focal_sites[0] - modified_left);
        ancestral_sequence.end = VariantIndex(last_focal_site + modified_right + 1);

        ancestral_sequence
    }

    /// Extend an ancestral sequence for a given set of sites (provided through an iterator).
    /// The ancestral state for each site is computed if the site is older than the focal sites for
    /// the inferred ancestral sequence.
    /// If `termination_condition` is true, the extension will be aborted once the set of samples
    /// believed to have derived their state from the ancestral sequence shrunk to half its original
    /// size. If it is false, the ancestral state for all sites will be calculated.
    ///
    /// # Parameter
    /// - `site_iter`: provides the sites for which to infer the common ancestor
    /// - `focal_site`: the site index on which the ancestor is based
    /// - `ancestral_sequence` the haplotype that is being generated
    /// - `termination_condition` if true, the extension will be terminated once enough samples
    /// diverge from the common focal site.
    ///
    /// # Returns
    /// Returns the number of continuous sites modified other than the focal site.
    fn extend_ancestor(
        &self,
        site_iter: &mut dyn Iterator<Item = (usize, &VariantSite)>,
        focal_site: usize,
        ancestral_sequence: &mut AncestralSequence,
        termination_condition: bool,
    ) -> usize {
        let mut modified_sites = 0;
        // the focal site is defined by its derived state
        ancestral_sequence[focal_site] = DERIVED_STATE;

        // the set of samples that are still considered part of the subtree derived from this ancestor
        // the generation process ends, once this set reaches half it's size
        let focal_set = self.sites[focal_site]
            .genotypes
            .iter()
            .enumerate()
            .filter(|&(_, s)| *s > 0)
            .map(|(i, _)| i)
            .enumerate()
            .collect::<Vec<_>>();
        let focal_set_size = focal_set.len();
        let focal_age = self.sites[focal_site].relative_age;

        // the current set of genotype states for the current site
        let mut current_set = vec![true; focal_set_size];

        // currently active samples. Initially, all samples that have the derived state at the focal site are active
        let mut active_samples_set = vec![true; focal_set_size];

        // marked for deletion from the active samples set. Deleted if a marked elements gets marked again
        let mut deletion_marks = vec![false; focal_set_size];

        // the size of the current set of samples
        let mut remaining_set_size = focal_set_size;

        for (variant_index, site) in site_iter {
            if site.relative_age > focal_age {
                if termination_condition {
                    // mask out the ones in the current site with the set of genotypes that derive from the ancestral sequence
                    let mut ones = 0;
                    focal_set.iter().for_each(|&(i, sample)| {
                        let state = site.genotypes[sample];
                        current_set[i] = state == DERIVED_STATE;

                        if active_samples_set[i] {
                            if state == DERIVED_STATE {
                                ones += 1;
                            }
                        }
                    });

                    // compute ancestral state
                    // tsinfer doesn't use integer division here, so subtraction is necessary
                    // for result parity
                    let consensus_state = if ones >= remaining_set_size - ones {
                        DERIVED_STATE
                    } else {
                        ANCESTRAL_STATE
                    };

                    modified_sites += 1;
                    ancestral_sequence[variant_index] = consensus_state;

                    // update the set of genotypes for next loop iteration
                    if consensus_state == DERIVED_STATE {
                        deletion_marks
                            .iter_mut()
                            .zip(current_set.iter_mut())
                            .for_each(|(mark, state)| {
                                *state = !*state;
                                *mark &= *state;
                            });
                    } else {
                        deletion_marks.iter_mut().zip(current_set.iter()).for_each(
                            |(mark, state)| {
                                *mark &= *state;
                            },
                        );
                    }

                    active_samples_set
                        .iter_mut()
                        .zip(deletion_marks.iter())
                        .for_each(|(active, mark)| {
                            *active &= !*mark;
                        });

                    remaining_set_size = active_samples_set.iter().filter(|&s| *s).count();

                    // we don't need the current set anymore, so we can reuse it for the next iteration
                    // as the deletion marks. We also don't need the deletion marks anymore, so we can
                    // overwrite them with the current set in the next iteration.
                    mem::swap(&mut deletion_marks, &mut current_set);

                    // tsinfer uses integer division for this test, so rounding errors are required
                    // for result parity
                    if remaining_set_size <= focal_set_size / 2 {
                        break;
                    }
                } else {
                    let mut ones = 0;
                    focal_set.iter().for_each(|&(_, sample)| {
                        if site.genotypes[sample] == DERIVED_STATE {
                            ones += 1;
                        }
                    });

                    // compute ancestral state
                    modified_sites += 1;
                    ancestral_sequence[variant_index] = if ones >= remaining_set_size - ones {
                        DERIVED_STATE
                    } else {
                        ANCESTRAL_STATE
                    };
                }
            } else {
                modified_sites += 1;
                ancestral_sequence[variant_index] = ANCESTRAL_STATE;
            }
        }

        modified_sites
    }

    /// Generate a set of ancestral sequences with the variant sites added to the generator.
    pub fn generate_ancestors(&self, sequence_length: SequencePosition) -> AncestorArray {
        // fixme this entire process is inefficient, we should sort the original sites
        let mut sites = self.sites.iter().enumerate().collect::<Vec<_>>();
        sites.sort_unstable_by(|(_, a), (_, b)| {
            a.relative_age.partial_cmp(&b.relative_age).unwrap()
        });

        let mut focal_sites: Vec<Vec<usize>> = Vec::new();
        let mut current_age: f64 = -1f64;

        // Todo we are reconstructing the hashmap a lot, this seems unnecessary
        let mut current_focal_sites: HashMap<Vec<u8>, Vec<usize>, BuildHasherDefault<XxHash64>> =
            Default::default();
        for (focal_site, site) in sites {
            if f64::abs(site.relative_age - current_age) < 1e-6 {
                if current_focal_sites.contains_key(&site.genotypes) {
                    current_focal_sites
                        .get_mut(&site.genotypes)
                        .expect("")
                        .push(focal_site);
                } else {
                    current_focal_sites.insert(site.genotypes.clone(), vec![focal_site]);
                }
            } else {
                current_age = site.relative_age;
                focal_sites.append(
                    &mut current_focal_sites
                        .iter()
                        .map(|(_, v)| v.clone())
                        .collect::<Vec<_>>(),
                );
                current_focal_sites = Default::default();
                current_focal_sites.insert(site.genotypes.clone(), vec![focal_site]);
            }
        }
        focal_sites.append(
            &mut current_focal_sites
                .iter()
                .map(|(_, v)| v.clone())
                .collect::<Vec<_>>(),
        );

        // break focal sites apart if they are interrupted by disagreeing old sites
        // (i.e. if they are from different subtrees in the ancestry)
        let mut broken_focal_sites = Vec::with_capacity((focal_sites.len() as f64 * 1.1) as usize);
        for mut focal_site in focal_sites {
            if focal_site.len() > 1 {
                focal_site.sort_unstable();
                let mut partial_focal_site = Vec::new();
                for i in 0..focal_site.len() - 1 {
                    partial_focal_site.push(focal_site[i]);
                    let focal_site_size = self.sites[focal_site[i]]
                        .genotypes
                        .iter()
                        .filter(|&&s| s == DERIVED_STATE)
                        .count();

                    let must_split = self.sites[focal_site[i] + 1..focal_site[i + 1]]
                        .iter()
                        .filter(|&site| site.relative_age > self.sites[focal_site[i]].relative_age)
                        .any(|site| {
                            let consensus = site
                                .genotypes
                                .iter()
                                .zip(self.sites[focal_site[i]].genotypes.iter())
                                .filter(|&(_, &focal_sample)| focal_sample == DERIVED_STATE)
                                .filter(|&(sample, _)| *sample == DERIVED_STATE)
                                .count();
                            consensus != focal_site_size && consensus != 0
                        });

                    if must_split {
                        broken_focal_sites.push(partial_focal_site);
                        partial_focal_site = Vec::new();
                    }
                }

                partial_focal_site.push(*focal_site.last().unwrap());
                broken_focal_sites.push(partial_focal_site);
            } else {
                broken_focal_sites.push(focal_site);
            }
        }

        let focal_sites = broken_focal_sites;

        // TODO make the parallelization optional (par_iter)
        let mut ancestors: Vec<_> = focal_sites
            .iter()
            .map(|focal_sites| self.generate_ancestor(focal_sites))
            .collect();

        // artificially add the root ancestor
        let mut ancestral_state = AncestralSequence::from_ancestral_state(self.sites.len(), 1.0);
        ancestral_state.end = VariantIndex(self.sites.len());
        ancestors.push(ancestral_state);

        // TODO parallel sort ancestors by age (par_sort_unstable_by)
        ancestors.sort_unstable_by(|a, b| {
            a.relative_age()
                .partial_cmp(&b.relative_age())
                .unwrap()
                .reverse()
        });

        AncestorArray::new(ancestors, self.variant_positions(), sequence_length)
    }

    /// Calculate the DNA sample sequences from the variant sites
    pub fn generate_samples(&self) -> Vec<AncestralSequence> {
        // transpose the variant sites matrix
        let mut samples = vec![
            AncestralSequence::from_ancestral_state(self.sites.len(), 0f64);
            self.sites[0].genotypes.len()
        ];
        for (site_index, genotypes) in self
            .sites
            .iter()
            .map(|variant_site| &variant_site.genotypes)
            .enumerate()
        {
            for (i, variant) in genotypes.iter().enumerate() {
                samples[i][site_index] = *variant;
            }
        }

        samples
    }

    pub fn tskit_export_sites(&self, path: &Path) -> io::Result<()> {
        let mut site_file = path.to_path_buf();
        site_file.push("sites.tsv");
        let mut writer = std::fs::File::create(site_file)?;
        writer.write_fmt(format_args!("id\tposition\tancestral_state\n"))?;

        for (i, site) in self.sites.iter().enumerate() {
            writer.write_fmt(format_args!(
                "{id}\t{position}\t{ancestral_state}\n",
                id = i,
                position = site.position,
                ancestral_state = 'C', // TODO encode actual ancestral state
            ))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ancestors::Ancestor;
    use crate::variants::VariantSite;

    #[test]
    fn compute_trivial_ancestors() {
        let site1 = vec![0, 0, 1, 0, 1];
        let site2 = vec![0, 1, 1, 0, 0];
        let site3 = vec![0, 1, 0, 0, 1];
        let site4 = vec![0, 0, 0, 1, 1];

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new_raw(site1, 1),
                VariantSite::new_raw(site2, 2),
                VariantSite::new_raw(site3, 3),
                VariantSite::new_raw(site4, 4),
            ],
        };

        let len = SequencePosition::from_usize(5);

        let ancestors = ag.generate_ancestors(len);

        assert_eq!(ancestors.len(), 5);

        // root ancestor
        assert_eq!(ancestors[Ancestor(0)].state, vec![0, 0, 0, 0]);

        assert!(ancestors.iter().any(|(_, a)| a.state == vec![1, 0, 0, 0]));
        assert!(ancestors.iter().any(|(_, a)| a.state == vec![0, 1, 0, 0]));
        assert!(ancestors.iter().any(|(_, a)| a.state == vec![0, 0, 1, 0]));
        assert!(ancestors.iter().any(|(_, a)| a.state == vec![0, 0, 0, 1]));
    }

    #[test]
    fn compute_multi_focal_ancestors() {
        let site1 = vec![0, 0, 0, 1, 1];
        let site2 = vec![0, 1, 1, 0, 0];
        let site3 = vec![0, 1, 1, 0, 0];
        let site4 = vec![0, 0, 0, 1, 1];

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new_raw(site1, 1),
                VariantSite::new_raw(site2, 2),
                VariantSite::new_raw(site3, 3),
                VariantSite::new_raw(site4, 4),
            ],
        };

        let len = SequencePosition::from_usize(5);

        let ancestors = ag.generate_ancestors(len);

        assert_eq!(ancestors.len(), 3);

        // root ancestor
        assert_eq!(ancestors[Ancestor(0)].state, vec![0, 0, 0, 0]);

        assert!(ancestors.iter().any(|(_, a)| a.state == vec![1, 0, 0, 1]));
        assert!(ancestors.iter().any(|(_, a)| a.state == vec![0, 1, 1, 0]));
    }
}
