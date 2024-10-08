use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;
use std::thread::available_parallelism;
use std::{io, mem};

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;
use twox_hash::XxHash64;

use crate::ancestors::{AncestorArray, AncestralSequence, ANCESTRAL_STATE, DERIVED_STATE};
use crate::variants::{VariantData, VariantIndex, VariantSequence, VariantSite};

/// Generates ancestral sequences for a given set of variant sites. The ancestral sequences are
/// generated using heuristic methods that use a small number of variant sites as focal sites and
/// infer the ancestral state for surrounding sites. For each set of focal sites, a single ancestral
/// sequence is generated.
pub struct AncestorGenerator {
    variant_data: VariantData,
    num_threads: u16,
}

impl AncestorGenerator {
    /// Create a new ancestor generator with given variant data. It will use the number of threads
    /// available on the system.
    pub fn from_variant_data(variant_data: VariantData) -> Self {
        let num_threads = available_parallelism()
            .unwrap_or_else(|err| {
                eprintln!(
                    "Error getting number of threads: {}. Defaulting to 1 thread.",
                    err
                );
                NonZeroUsize::new(1).unwrap()
            })
            .get() as u16;

        Self::with_parallelism(variant_data, num_threads)
    }

    /// Create a new ancestor generator with given variant data and number of threads to use.
    pub fn with_parallelism(variant_data: VariantData, num_threads: u16) -> Self {
        Self {
            variant_data,
            num_threads,
        }
    }

    /// For a given set of focal sites, compute an ancestor that uses those focal sites. The focal
    /// site is given by a sorted set of indices into the set of variant sites.
    fn generate_ancestor(
        &self,
        buffer_sequence: &mut VariantSequence,
        focal_sites: &[VariantIndex],
    ) -> AncestralSequence {
        debug_assert!(!focal_sites.is_empty());
        debug_assert!(focal_sites.windows(2).all(|sites| sites[0] < sites[1]));

        let sites = self.variant_data.len();

        // extend ancestor to the left of the first focal site
        let modified_left = self.extend_ancestor(
            &mut self
                .variant_data
                .iter()
                .enumerate()
                .rev()
                .skip(sites - focal_sites[0].unwrap() as usize),
            focal_sites[0],
            buffer_sequence,
            true,
        );

        // infer ancestor between focal sites.
        for foc_index in 0..focal_sites.len() - 1 {
            let focal_site_i = focal_sites[foc_index];
            let focal_site_j = focal_sites[foc_index + 1];
            self.extend_ancestor(
                &mut self
                    .variant_data
                    .iter()
                    .enumerate()
                    .skip((focal_site_i.unwrap() + 1) as usize)
                    .take((focal_site_j.unwrap() - focal_site_i.unwrap() - 1) as usize),
                focal_site_i,
                buffer_sequence,
                false,
            );
        }

        // extend ancestor to the right of the last focal site
        let last_focal_site = *focal_sites.last().unwrap();
        let modified_right = self.extend_ancestor(
            &mut self
                .variant_data
                .iter()
                .enumerate()
                .skip((last_focal_site.unwrap() + 1) as usize),
            last_focal_site,
            buffer_sequence,
            true,
        );

        let ancestral_sequence = AncestralSequence::copy_from(
            buffer_sequence,
            focal_sites[0] - modified_left,
            last_focal_site + modified_right + 1,
            focal_sites.iter().copied().collect(),
            self.variant_data[focal_sites[0]].relative_age,
        );

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
        focal_site: VariantIndex,
        ancestral_sequence: &mut VariantSequence,
        termination_condition: bool,
    ) -> u32 {
        let mut modified_sites = 0;
        // the focal site is defined by its derived state
        ancestral_sequence[focal_site] = DERIVED_STATE;

        // the set of samples that are still considered part of the subtree derived from this ancestor
        // the generation process ends, once this set reaches half it's size
        let focal_set = self.variant_data[focal_site]
            .genotypes
            .iter()
            .enumerate()
            .filter(|&(_, s)| *s > 0)
            .map(|(i, _)| i)
            .enumerate()
            .collect::<Vec<_>>();
        let focal_set_size = focal_set.len();
        let focal_age = self.variant_data[focal_site].relative_age;

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
    pub fn generate_ancestors(self) -> AncestorArray {
        // we need to clone and sort the sites by age, which we cannot do for the original sites,
        // because they are sorted by sequence position, which we need to preserve
        let mut sites = self.variant_data.iter_with_index().collect::<Vec<_>>();

        if self.num_threads == 1 {
            sites.sort_unstable_by(|(_, a), (_, b)| {
                a.relative_age.partial_cmp(&b.relative_age).unwrap()
            });
        } else {
            // we assume the global threadpool is already initialized with the correct thread count
            sites.par_sort_unstable_by(|(_, a), (_, b)| {
                a.relative_age.partial_cmp(&b.relative_age).unwrap()
            });
        }

        let mut focal_sites: Vec<Vec<VariantIndex>> = Vec::new();
        let mut current_age: f64 = -1f64;

        // Todo we are reconstructing the hashmap a lot, this seems unnecessary
        // Todo parallelize this
        let mut current_focal_sites: HashMap<
            Vec<u8>,
            Vec<VariantIndex>,
            BuildHasherDefault<XxHash64>,
        > = Default::default();
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
                    let focal_site_size = self.variant_data[focal_site[i]]
                        .genotypes
                        .iter()
                        .filter(|&&s| s == DERIVED_STATE)
                        .count();

                    let must_split = self.variant_data[focal_site[i] + 1..focal_site[i + 1]]
                        .iter()
                        .filter(|&site| {
                            site.relative_age > self.variant_data[focal_site[i]].relative_age
                        })
                        .any(|site| {
                            let consensus = site
                                .genotypes
                                .iter()
                                .zip(self.variant_data[focal_site[i]].genotypes.iter())
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

        let mut ancestors: Vec<_> = if self.num_threads == 1 {
            let mut buffer = VariantSequence::from_ancestral_state(self.variant_data.len() as u32);
            focal_sites
                .iter()
                .map(|focal_sites| self.generate_ancestor(&mut buffer, focal_sites))
                .collect()
        } else {
            // we assume the global threadpool is already initialized with the correct thread count
            focal_sites
                .par_iter()
                .map_with(
                    VariantSequence::from_ancestral_state(self.variant_data.len() as u32),
                    |buffer, focal_sites| self.generate_ancestor(buffer, focal_sites),
                )
                .collect()
        };

        // artificially add the root ancestor
        let mut ancestral_state =
            AncestralSequence::from_ancestral_state(self.variant_data.len() as u32, 2.0);
        ancestral_state.end = VariantIndex(self.variant_data.len() as u32);
        ancestors.push(ancestral_state);

        if self.num_threads == 1 {
            ancestors.sort_unstable_by(|a, b| {
                a.relative_age()
                    .partial_cmp(&b.relative_age())
                    .unwrap()
                    .reverse()
            });
        } else {
            ancestors.par_sort_unstable_by(|a, b| {
                a.relative_age()
                    .partial_cmp(&b.relative_age())
                    .unwrap()
                    .reverse()
            });
        }

        AncestorArray::new(ancestors, self.variant_data)
    }

    /// Calculate the DNA sample sequences from the variant sites
    pub fn generate_samples(&self) -> Vec<AncestralSequence> {
        // transpose the variant sites matrix
        let mut samples =
            vec![
                AncestralSequence::from_ancestral_state(self.variant_data.len() as u32, 0f64);
                self.variant_data.get_num_samples() as usize
            ];
        for (site_index, genotypes) in self
            .variant_data
            .iter()
            .map(|variant_site| &variant_site.genotypes)
            .enumerate()
        {
            for (i, variant) in genotypes.iter().enumerate() {
                samples[i][VariantIndex::from_usize(site_index)] = *variant;
            }
        }

        samples
    }

    pub fn tskit_export_sites(&self, path: &Path) -> io::Result<()> {
        let mut site_file = path.to_path_buf();
        site_file.push("sites.tsv");
        let mut writer = std::fs::File::create(site_file)?;
        writer.write_fmt(format_args!("id\tposition\tancestral_state\n"))?;

        for (i, site) in self.variant_data.iter().enumerate() {
            writer.write_fmt(format_args!(
                "{id}\t{position}\t{ancestral_state}\n",
                id = i,
                position = site.position,
                ancestral_state = site.ancestral_state,
            ))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::Ancestor;
    use crate::variants::VariantDataBuilder;

    use super::*;

    #[test]
    fn compute_trivial_ancestors() {
        let site1 = vec![0, 0, 1, 0, 1];
        let site2 = vec![0, 1, 1, 0, 0];
        let site3 = vec![0, 1, 0, 0, 1];
        let site4 = vec![0, 0, 0, 1, 1];

        let variant_data = VariantDataBuilder::from_iter(
            5,
            vec![
                (site1, 1, 'G', 'A'),
                (site2, 2, 'G', 'A'),
                (site3, 3, 'G', 'C'),
                (site4, 4, 'G', 'T'),
            ],
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();

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

        let variant_data = VariantDataBuilder::from_iter(
            5,
            vec![
                (site1, 1, 'G', 'A'),
                (site2, 2, 'G', 'C'),
                (site3, 3, 'G', 'T'),
                (site4, 4, 'G', 'T'),
            ],
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();

        assert_eq!(ancestors.len(), 3);

        // root ancestor
        assert_eq!(ancestors[Ancestor(0)].state, vec![0, 0, 0, 0]);

        assert!(ancestors.iter().any(|(_, a)| a.state == vec![1, 0, 0, 1]));
        assert!(ancestors.iter().any(|(_, a)| a.state == vec![0, 1, 1, 0]));
    }
}
