use std::collections::HashMap;
use crate::dna::VariantSite;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;
use std::fmt::{Debug, Formatter};
use std::mem;

const ANCESTRAL_STATE: u8 = 0;
const DERIVED_STATE: u8 = 1;

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
// TODO: is this enough? If it isn't, this must be replaced with a packed sequence coding for whatever
//  state we need, ideally with a optimization because we expect most of the entries to reference the
//  ancestral state
#[derive(Clone)]
pub struct AncestralSequence {
    pub(crate) state: Vec<u8>,
    pub(crate) focal_sites: Vec<usize>,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) age: f64,
}

impl AncestralSequence {
    // TODO shouldnt be public
    pub fn from_ancestral_state(len: usize, age: f64) -> Self {
        AncestralSequence {
            state: vec![0u8; len],
            focal_sites: Vec::new(),
            start: 0,
            end: 0,
            age,
        }
    }

    fn set_unchecked(&mut self, index: usize, value: u8) {
        self.state[index] = value;
    }

    pub fn len(&self) -> usize {
        self.state.len()
    }
}

impl Debug for AncestralSequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("AncestralSequence { ")?;
        f.write_fmt(format_args!("focal_sites={:?},\t", self.focal_sites))?;

        f.write_str("genotype=[ ")?;
        let mut iter = self.state.iter();
        let mut idx = 0;
        while let Some(b) = iter.next() {
            if idx < self.start || idx > self.end {
                f.write_str("-")?;
            } else {
                f.write_fmt(format_args!("{}", b))?;
            }

            if iter.len() > 0 {
                f.write_str(", ")?;
            }

            idx += 1;
        }
        f.write_str(" ] }")?;
        Ok(())
    }
}

pub struct AncestorGenerator {
    sites: Vec<VariantSite>,
}

impl AncestorGenerator {
    pub fn new(sites: Vec<VariantSite>) -> Self {
        Self {
            sites
        }
    }

    /// For a given set of focal sites, compute an ancestor that uses those focal sites.
    pub fn generate_ancestor(&self, focal_sites: &[usize]) -> AncestralSequence {
        debug_assert!(!focal_sites.is_empty());
        debug_assert!(focal_sites.windows(2).all(|sites| sites[0] < sites[1]));

        let mut ancestral_sequence = AncestralSequence::from_ancestral_state(
            self.sites.len(),
            self.sites[focal_sites[0]].relative_age,
        );
        let sites = self.sites.len();

        // extend ancestor to the left of the first focal site
        let start = focal_sites[0]
            - self.extend_ancestor(
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
        for focal_site in 1..focal_sites.len() - 1 {
            let focal_site_i = focal_sites[focal_site];
            let focal_site_j = focal_sites[focal_site + 1];
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
        let end = last_focal_site
            + self.extend_ancestor(
                &mut self.sites.iter().enumerate().skip(last_focal_site + 1),
                last_focal_site,
                &mut ancestral_sequence,
                true,
            );

        // TODO truncate ancestral sequence to save memory. this should be implemented using functionality from BitVec
        ancestral_sequence.focal_sites = focal_sites.to_vec();
        ancestral_sequence.start = start;
        ancestral_sequence.end = end;

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
        ancestral_sequence.set_unchecked(focal_site, DERIVED_STATE);

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
                    let consensus_state = if ones >= remaining_set_size - ones {
                        DERIVED_STATE
                    } else {
                        ANCESTRAL_STATE
                    };

                    modified_sites += 1;
                    ancestral_sequence.set_unchecked(variant_index, consensus_state);

                    // update the set of genotypes for next loop iteration
                    if consensus_state == DERIVED_STATE {
                        deletion_marks.iter_mut().zip(current_set.iter_mut()).for_each(
                            |(mark, state)| {
                                *state = !*state;
                                *mark &= *state;
                            },
                        );
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

                    if remaining_set_size <= focal_set_size / 2 {
                        break;
                    }
                } else {
                    let mut ones = 0;
                    focal_set.iter().for_each(|&(i, sample)| {
                        if site.genotypes[sample] == DERIVED_STATE {
                            ones += 1;
                        }
                    });

                    // compute ancestral state
                    modified_sites += 1;
                    if ones >= remaining_set_size - ones {
                        ancestral_sequence.set_unchecked(variant_index, DERIVED_STATE);
                    } else {
                        ancestral_sequence.set_unchecked(variant_index, ANCESTRAL_STATE);
                    };
                }
            } else {
                modified_sites += 1;
                ancestral_sequence.set_unchecked(variant_index, ANCESTRAL_STATE);
            }
        }

        modified_sites
    }

    pub fn generate_ancestors(&self) -> Vec<AncestralSequence> {
        // TODO we have to sort focal sites by time and group those with equal genotype distributions
        //  together

        // fixme this entire process is inefficient, we should sort the original sites
        let mut sites = self.sites.iter().enumerate().collect::<Vec<_>>();
        sites.sort_unstable_by(|(_, a), (_, b)| a.relative_age.partial_cmp(&b.relative_age).unwrap());

        let mut focal_sites: Vec<Vec<usize>> = Vec::new();
        let mut current_age: f64 = -1f64;
        let mut current_focal_sites = HashMap::<Vec<u8>, Vec<usize>>::new();
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
                current_focal_sites = HashMap::new();
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

        // TODO make the parallelization optional

        let ancestors = focal_sites
            .par_iter()
            .map(|focal_sites| self.generate_ancestor(focal_sites))
            .collect();

        // TODO sort ancestors by age

        ancestors
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::{AncestorGenerator, DERIVED_STATE};
    use crate::dna::VariantSite;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::time::Instant;

    #[test]
    fn compute_trivial_ancestors() {
        let site1 = vec![0, 0, 1, 0, 1];
        let site2 = vec![0, 1, 1, 0, 0];
        let site3 = vec![0, 1, 0, 0, 1];
        let site4 = vec![0, 0, 0, 1, 1];

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site3, 3),
                VariantSite::new(site4, 4),
            ],
        };

        let ancestors = ag.generate_ancestors();

        // TODO find out why tsinfer generates the root ancestor twice
        // TODO generate root ancestor
        assert_eq!(ancestors.len(), 4); // TODO 5

        assert!(ancestors.iter().any(|a| a.state == vec![1, 0, 0, 0]));
        assert!(ancestors.iter().any(|a| a.state == vec![0, 1, 0, 0]));
        assert!(ancestors.iter().any(|a| a.state == vec![0, 0, 1, 0]));
        assert!(ancestors.iter().any(|a| a.state == vec![0, 0, 0, 1]));
    }

    #[test]
    fn compute_multi_focal_ancestors() {
        let site1 = vec![0, 0, 0, 1, 1];
        let site2 = vec![0, 1, 1, 0, 0];
        let site3 = vec![0, 1, 1, 0, 0];
        let site4 = vec![0, 0, 0, 1, 1];

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site3, 3),
                VariantSite::new(site4, 4),
            ],
        };

        let ancestors = ag.generate_ancestors();

        // TODO find out why tsinfer generates the root ancestor twice
        // TODO generate root ancestor
        assert_eq!(ancestors.len(), 2); // TODO 3

        assert!(ancestors.iter().any(|a| a.state == vec![1, 0, 0, 1]));
        assert!(ancestors.iter().any(|a| a.state == vec![0, 1, 1, 0]));
    }

    /// Reads a specially formatted text file that contains data about variant sites intended for
    /// unit testing. The data was generated by dumping it from tsinfer.
    fn read_variant_dump(path: &str) -> Vec<VariantSite> {
        let input = File::open(path).expect("could not find test data");
        let reader = BufReader::new(input);
        reader
            .lines()
            .enumerate()
            .map(|(pos, line)| {
                VariantSite::new(
                    line.expect("unexpected io error")
                        .trim()
                        .split(" ")
                        .map(|s| s.parse().expect("corrupt input data"))
                        .collect::<Vec<_>>(),
                    pos,
                )
            })
            // filter out singleton sites TODO: this should be done by a builder interface
            .filter(|seq| {
                seq.genotypes
                    .iter()
                    .filter(|&state| *state == DERIVED_STATE)
                    .count()
                    > 1
            })
            .collect::<Vec<_>>()
    }

    /// Reads a specially formatted text file that contains data about ancestral sequences intended
    /// for unit testing. The data was generated by dumping it from tsinfer.
    fn read_ancestor_dump(path: &str, sequences: usize, sequence_length: usize) -> Vec<Vec<u8>> {
        let input = File::open(path).expect("could not find test results");
        let reader = BufReader::new(input);
        let mut tsinfer_ancestors = vec![vec![0; sequence_length]; sequences];
        reader.lines().for_each(|line| {
            let line = line.expect("unexpected IO error");
            let mut line_parts = line.splitn(2, ": ");
            let index: usize = line_parts
                .next()
                .expect("corrupted test results")
                .parse()
                .expect("corrupted test results");
            line_parts
                .next()
                .expect("corrupted test results")
                .trim()
                .split(" ")
                .enumerate()
                .for_each(|(j, bit)| {
                    tsinfer_ancestors[index][j] = bit.parse().expect("corrupted test results")
                });
        });
        tsinfer_ancestors
    }

    #[test]
    fn compute_chr20_40_variants() {
        let variant_sites = read_variant_dump("testdata/chr20_40variants.txt");
        // according to tsinfer, we have 22 variant sites
        assert_eq!(variant_sites.len(), 22);

        let ag = AncestorGenerator {
            sites: variant_sites,
        };

        let ancestors = ag.generate_ancestors();

        let tsinfer_ancestors = read_ancestor_dump("testdata/chr20_40ancestors.txt", 22, 22);

        for (index, ancestor) in ancestors.iter().enumerate() {
            assert_eq!(ancestor.focal_sites.len(), 1);
            for (pos, &state) in ancestor.state.iter().enumerate() {
                assert_eq!(state, tsinfer_ancestors[ancestor.focal_sites[0]][pos]);
            }
        }
    }

    #[test]
    fn compute_chr20_10k_variants() {
        let variant_sites = read_variant_dump("testdata/chr20_10k_variants.txt");

        let ag = AncestorGenerator {
            sites: variant_sites,
        };

        let start = Instant::now();
        let ancestors = ag.generate_ancestors();
        println!("time passed: {:?}", start.elapsed());

        // TODO check with tsinfer results
    }
}
