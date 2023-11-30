use crate::dna::AncestralSequence;
use crate::samples::VariantSite;

const ANCESTRAL_STATE: u64 = 0;
const DERIVED_STATE: u64 = 1;

pub struct AncestorGenerator {
    sites: Vec<VariantSite>,
}

impl AncestorGenerator {
    fn generate_ancestor(&self, focal_site: usize) -> AncestralSequence {
        let mut ancestral_sequence = AncestralSequence::from_ancestral_state(self.sites.len());
        let sites = self.sites.len();
        ancestral_sequence.set_unchecked(focal_site, 1);

        self.extend_ancestor(&mut self.sites.iter().enumerate().skip(focal_site + 1), focal_site, &mut ancestral_sequence);
        self.extend_ancestor(&mut self.sites.iter().enumerate().rev().skip(sites - focal_site), focal_site, &mut ancestral_sequence);

        // TODO truncate ancestral sequence to save memory. this should be implemented using functionality from BitVec

        ancestral_sequence
    }

    fn extend_ancestor(&self, site_iter: &mut dyn Iterator<Item=(usize, &VariantSite)>, focal_site: usize, ancestral_sequence: &mut AncestralSequence) {
        // the set of samples that are still considered part of the subtree derived from this ancestor
        // the generation process ends, once this set reaches half it's size
        let mut derived_set = self.sites[focal_site].genotypes.clone();
        let focal_site_age = self.sites[focal_site].relative_age;
        let derived_set_size = derived_set.count_ones();

        // the size of the current set of samples
        let mut current_set_size = derived_set_size;

        for (variant_index, site) in site_iter {
            if site.relative_age >= focal_site_age {
                let masked_set = site.genotypes.mask_and(&derived_set).expect("didn't expect different length variant sites").to_bit_vec();
                let ancestral_state = if masked_set.count_ones() > current_set_size / 2 { 1 } else { 0 };
                ancestral_sequence.set_unchecked(variant_index, ancestral_state);
                derived_set = masked_set;
                current_set_size = derived_set.count_ones();
                // todo make sure the proper termination condition is used
                if current_set_size < derived_set_size / 2 {
                    // todo remember where we stopped the ancestor
                    break;
                }
            }
        }
    }

    pub fn generate_ancestors(&self) -> Vec<AncestralSequence> {
        let mut ancestors = Vec::with_capacity(self.sites.len());

        for (focal_site, _) in self.sites.iter().enumerate() {
            let ancestral_sequence = self.generate_ancestor(focal_site);
            println!("Focal Site {}: {:?}", focal_site, ancestral_sequence);
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
                    BitVec::from_bits(&line.expect("unexpected io error")
                        .trim()
                        .split(" ")
                        .map(|s| s.parse().expect("corrupt input data"))
                        .collect::<Vec<_>>()),
                    pos,
                )
            })
            // filter out singleton sites TODO: this should be done by a builder interface
            .filter(|seq| seq.genotypes.iter().filter(|&state| state == DERIVED_STATE).count() > 1)
            .collect::<Vec<_>>();

        // according to tsinfer, we have 22 variant sites
        assert_eq!(variant_sites.len(), 22);

        let ag = AncestorGenerator {
            sites: variant_sites
        };

        let ancestors = ag.generate_ancestors();
    }
}
