use std::fmt::{Debug, Formatter};
use crate::samples::VariantSite;
use std::mem;
use vers_vecs::BitVec;

const ANCESTRAL_STATE: u64 = 0;
const DERIVED_STATE: u64 = 1;

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
// TODO: is this enough? If it isn't, this must be replaced with a packed sequence coding for whatever
//  state we need, ideally with a optimization because we expect most of the entries to reference the
//  ancestral state
pub struct AncestralSequence {
    state: BitVec,
    start: usize,
    end: usize,
    age: f64,
}

impl AncestralSequence {
    fn from_ancestral_state(len: usize, age: f64) -> Self {
        AncestralSequence {
            state: BitVec::from_zeros(len),
            start: 0,
            end: 0,
            age,
        }
    }

    fn set_unchecked(&mut self, index: usize, value: u64) {
        self.state.set_unchecked(index, value);
    }
}

impl Debug for AncestralSequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("AncestralSequence [")?;
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
        f.write_str("]")?;
        Ok(())
    }
}

pub struct AncestorGenerator {
    sites: Vec<VariantSite>,
}

impl AncestorGenerator {
    /// For a given set of focal sites, compute an ancestor that uses those focal sites.
    pub fn generate_ancestor(&self, focal_sites: &[usize]) -> AncestralSequence {
        debug_assert!(!focal_sites.is_empty());
        debug_assert!(focal_sites.windows(2).all(|sites| sites[0] < sites[1]));

        let mut ancestral_sequence = AncestralSequence::from_ancestral_state(self.sites.len(), self.sites[focal_sites[0]].relative_age);
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
        let focal_set = self.sites[focal_site].genotypes.iter1().enumerate().collect::<Vec<_>>();
        let focal_set_size = self.sites[focal_site].genotypes.rank1(self.sites[focal_site].genotypes.len());
        let focal_age = self.sites[focal_site].relative_age;

        // the current set of genotype states for the current site
        let mut current_set = BitVec::from_ones(focal_set_size);

        // currently active samples. Initially, all samples that have the derived state at the focal site are active
        let mut active_samples_set = BitVec::from_ones(focal_set_size);

        // marked for deletion from the active samples set. Deleted if a marked elements gets marked again
        let mut deletion_marks = BitVec::from_zeros(focal_set_size);

        // the size of the current set of samples
        let mut remaining_set_size = focal_set_size;

        for (variant_index, site) in site_iter {
            if site.relative_age > focal_age {
                if termination_condition {
                    // mask out the ones in the current site with the set of genotypes that derive from the ancestral sequence
                    let mut ones = 0;
                    focal_set.iter().for_each(|&(i, sample)| {
                        let state = site.genotypes.get_unchecked(sample);
                        current_set.set_unchecked(i, state);

                        if active_samples_set.get_unchecked(i) == 1 {
                            if state == DERIVED_STATE {
                                ones += 1;
                            }
                        }
                    });

                    // compute ancestral state
                    let consensus_state = if ones >= remaining_set_size / 2 {
                        DERIVED_STATE
                    } else {
                        ANCESTRAL_STATE
                    };

                    modified_sites += 1;
                    ancestral_sequence.set_unchecked(variant_index, consensus_state);

                    // update the set of genotypes for next loop iteration
                    if consensus_state == DERIVED_STATE {
                        deletion_marks
                            .apply_mask_custom(&current_set, |a, b| a & !b)
                            .expect("expect same length variant sites")
                    } else {
                        deletion_marks
                            .apply_mask_custom(&current_set, |a, b| a & b)
                            .expect("expect same length variant sites")
                    }
                    active_samples_set
                        .apply_mask_custom(&deletion_marks, |a, b| a & !b)
                        .expect("expect same length variant sites");

                    remaining_set_size = active_samples_set.count_ones() as usize;

                    // we don't need the current set anymore, so we can reuse it for the next iteration
                    // as the deletion marks. We also don't need the deletion marks anymore, so we can
                    // overwrite them with the current set in the next iteration.
                    mem::swap(&mut deletion_marks, &mut current_set);

                    if remaining_set_size < focal_set_size / 2 {
                        break;
                    }
                } else {
                    // TODO implement using the new iterators instead of masking the whole set
                    // let masked_set = site
                    //     .genotypes
                    //     .mask_and(&active_set)
                    //     .expect("didn't expect different length variant sites");
                    // let ancestral_state = if masked_set.count_ones() > remaining_set_size / 2 {
                    //     1
                    // } else {
                    //     0
                    // };

                    // modified_sites += 1;
                    // ancestral_sequence.set_unchecked(variant_index, ancestral_state);
                }
            } else {
                modified_sites += 1;
                ancestral_sequence.set_unchecked(variant_index, ANCESTRAL_STATE);
            }
        }

        modified_sites
    }

    pub fn generate_ancestors(&self) -> Vec<AncestralSequence> {
        let mut ancestors = Vec::with_capacity(self.sites.len());

        // TODO we have to sort focal sites by time and group those with equal genotype distributions
        //  together

        // TODO we have to break apart ancestors that have older sites in between them where not all
        //  samples derived from the ancestor agree on the state (I'm unsure why though)

        for (focal_site, _) in self.sites.iter().enumerate() {
            let ancestral_sequence = self.generate_ancestor(&[focal_site]);
            // println!(
            //     "Focal Site {} (time: {}): {:?}",
            //     focal_site, self.sites[focal_site].relative_age, ancestral_sequence
            // );
            ancestors.push(ancestral_sequence);
        }

        ancestors
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::{AncestorGenerator, DERIVED_STATE};
    use crate::samples::VariantSite;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::time::Instant;
    use vers_vecs::BitVec;

    #[test]
    fn compute_trivial_ancestors() {
        let site1 = BitVec::from_bits(&[0, 0, 1, 0, 1]);
        let site2 = BitVec::from_bits(&[0, 1, 1, 0, 0]);
        let site3 = BitVec::from_bits(&[0, 1, 0, 0, 1]);
        let site4 = BitVec::from_bits(&[0, 0, 0, 1, 1]);

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

        assert_eq!(ancestors[0].state, BitVec::from_bits(&[1, 0, 0, 0]));
        assert_eq!(ancestors[1].state, BitVec::from_bits(&[0, 1, 0, 0]));
        assert_eq!(ancestors[2].state, BitVec::from_bits(&[0, 0, 1, 0]));
        assert_eq!(ancestors[3].state, BitVec::from_bits(&[0, 0, 0, 1]));
    }

    #[test]
    fn compute_multi_focal_ancestors() {
        let site1 = BitVec::from_bits(&[0, 0, 0, 1, 1]);
        let site2 = BitVec::from_bits(&[0, 1, 1, 0, 0]);
        let site3 = BitVec::from_bits(&[0, 1, 1, 0, 0]);
        let site4 = BitVec::from_bits(&[0, 0, 0, 1, 1]);

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

        assert_eq!(ancestors[1].state, BitVec::from_bits(&[0, 1, 1, 0]));
        assert_eq!(ancestors[2].state, BitVec::from_bits(&[1, 0, 0, 1]));
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
                    BitVec::from_bits(
                        &line
                            .expect("unexpected io error")
                            .trim()
                            .split(" ")
                            .map(|s| s.parse().expect("corrupt input data"))
                            .collect::<Vec<_>>(),
                    ),
                    pos,
                )
            })
            // filter out singleton sites TODO: this should be done by a builder interface
            .filter(|seq| {
                seq.genotypes
                    .iter()
                    .filter(|&state| state == DERIVED_STATE)
                    .count()
                    > 1
            })
            .collect::<Vec<_>>()
    }

    /// Reads a specially formatted text file that contains data about ancestral sequences intended
    /// for unit testing. The data was generated by dumping it from tsinfer.
    fn read_ancestor_dump(path: &str, sequences: usize, sequence_length: usize) -> Vec<BitVec> {
        let input = File::open(path).expect("could not find test results");
        let reader = BufReader::new(input);
        let mut tsinfer_ancestors = vec![BitVec::from_zeros(sequence_length); sequences];
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
                    tsinfer_ancestors[index]
                        .set_unchecked(j, bit.parse().expect("corrupted test results"))
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
            for (pos, state) in ancestor.state.iter().enumerate() {
                assert_eq!(state, tsinfer_ancestors[index].get_unchecked(pos));
            }
        }
    }

    // #[test]
    fn compute_chr20_10k_variants() {
        let variant_sites = read_variant_dump("testdata/chr20_10k_variants.txt");

        let ag = AncestorGenerator {
            sites: variant_sites,
        };

        let start = Instant::now();
        let ancestors = ag.generate_ancestors();
        println!("time passed: {:?}", start.elapsed());

        // let tsinfer_ancestors = read_ancestor_dump("testdata/chr20_10k_ancestors.txt", 5177, -1);
        //
        // for (index, ancestor) in ancestors.iter().enumerate() {
        //     for (pos, state) in ancestor.state.iter().enumerate() {
        //         assert_eq!(state, tsinfer_ancestors[index].get_unchecked(pos));
        //     }
        // }
    }
}
