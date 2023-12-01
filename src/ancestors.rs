use crate::dna::AncestralSequence;
use crate::samples::VariantSite;
use vers_vecs::BitVec;

const ANCESTRAL_STATE: u64 = 0;
const DERIVED_STATE: u64 = 1;

pub struct AncestorGenerator {
    sites: Vec<VariantSite>,
}

impl AncestorGenerator {
    /// For a given set of focal sites, compute an ancestor that uses those focal sites.
    fn generate_ancestor(&self, focal_sites: &[usize]) -> AncestralSequence {
        debug_assert!(!focal_sites.is_empty());
        debug_assert!(focal_sites.windows(2).all(|sites| sites[0] < sites[1]));

        let mut ancestral_sequence = AncestralSequence::from_ancestral_state(self.sites.len());
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
        site_iter: &mut dyn Iterator<Item=(usize, &VariantSite)>,
        focal_site: usize,
        ancestral_sequence: &mut AncestralSequence,
        termination_condition: bool,
    ) -> usize {
        let mut modified_sites = 0;
        // the focal site is defined by its derived state
        ancestral_sequence.set_unchecked(focal_site, DERIVED_STATE);

        // the set of samples that are still considered part of the subtree derived from this ancestor
        // the generation process ends, once this set reaches half it's size
        let mut active_set = self.sites[focal_site].genotypes.clone();
        let mut previous_deletion_marks = BitVec::from_zeros(active_set.len());

        let focal_site_age = self.sites[focal_site].relative_age;
        let active_set_start_size = active_set.count_ones();

        // the size of the current set of samples
        let mut remaining_set_size = active_set_start_size;

        for (variant_index, site) in site_iter {
            if site.relative_age > focal_site_age {
                if termination_condition {
                    // mask out the ones in the current site with the set of genotypes that derive from the ancestral sequence
                    let masked_set = site
                        .genotypes
                        .mask_and(&active_set)
                        .expect("expect same length variant sites")
                        .to_bit_vec();

                    // compute ancestral state
                    let masked_ones = masked_set.count_ones();
                    let ancestral_state = if masked_ones >= remaining_set_size / 2 {
                        DERIVED_STATE
                    } else {
                        ANCESTRAL_STATE
                    };

                    modified_sites += 1;
                    ancestral_sequence.set_unchecked(variant_index, ancestral_state);

                    // update the set of genotypes for next loop iteration
                    // TODO implement in-place masking in BitVec
                    let deletion_marks = if ancestral_state == DERIVED_STATE {
                        active_set
                            .mask_custom(&masked_set, |a, b| (a ^ b) & a)
                            .expect("expect same length variant sites")
                            .to_bit_vec()
                    } else {
                        active_set
                            .mask_custom(&masked_set, |a, b| a & b)
                            .expect("expect same length variant sites")
                            .to_bit_vec()
                    };

                    let delete_set = deletion_marks
                        .mask_and(&previous_deletion_marks)
                        .expect("expect same length variant sites")
                        .to_bit_vec();

                    // TODO modify in-place
                    let new_active_set = active_set
                        .mask_custom(&delete_set, |a, b| a & !b)
                        .expect("expect same length variant sites")
                        .to_bit_vec();

                    active_set = new_active_set;

                    remaining_set_size = active_set.count_ones();

                    previous_deletion_marks = deletion_marks;

                    if remaining_set_size < active_set_start_size / 2 {
                        break;
                    }
                } else {
                    let masked_set = site
                        .genotypes
                        .mask_and(&active_set)
                        .expect("didn't expect different length variant sites");
                    let ancestral_state = if masked_set.count_ones() > remaining_set_size / 2 {
                        1
                    } else {
                        0
                    };

                    modified_sites += 1;
                    ancestral_sequence.set_unchecked(variant_index, ancestral_state);
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
            println!(
                "Focal Site {} (time: {}): {:?}",
                focal_site, self.sites[focal_site].relative_age, ancestral_sequence
            );
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
    use vers_vecs::BitVec;

    #[test]
    fn compute_trivial_ancestors() {
        let site1 = BitVec::from_bits(&[0, 0, 1, 0, 1, 0]);
        let site2 = BitVec::from_bits(&[0, 1, 1, 0, 0, 0]);
        let site3 = BitVec::from_bits(&[0, 1, 0, 0, 1, 0]);
        let site4 = BitVec::from_bits(&[0, 0, 0, 1, 1, 0]);

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

    #[test]
    fn compute_chr20_40_variants() {
        let input = File::open("testdata/chr20_40variants.txt").expect("could not find test data");
        let reader = BufReader::new(input);
        let variant_sites = reader
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
            .collect::<Vec<_>>();

        // according to tsinfer, we have 22 variant sites
        assert_eq!(variant_sites.len(), 22);

        let ag = AncestorGenerator {
            sites: variant_sites,
        };

        let ancestors = ag.generate_ancestors();
    }
}
