use crate::dna::AncestralSequence;
use crate::samples::VariantSite;

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
    use vers_vecs::BitVec;
    use crate::ancestors::AncestorGenerator;
    use crate::samples::VariantSite;

    #[test]
    fn compute_trivial_ancestors() {
        let gen2 = BitVec::from_bits(&[0, 0, 1, 0, 1, 0]);
        let gen3 = BitVec::from_bits(&[0, 1, 1, 0, 0, 0]);
        let gen4 = BitVec::from_bits(&[0, 1, 0, 0, 1, 0]);
        let gen5 = BitVec::from_bits(&[0, 0, 0, 1, 1, 0]);

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new(gen2, 1),
                VariantSite::new(gen3, 2),
                VariantSite::new(gen4, 3),
                VariantSite::new(gen5, 4),
            ],
        };

        let ancestors = ag.generate_ancestors();

        // TODO find out why tsinfer generates the root ancestor twice
        // TODO generate root ancestor
        assert_eq!(ancestors.len(), 4); // TODO 5

        // TODO bitvec needs Eq implementation
        // assert_eq!(ancestors[0].state, BitVec::from_bits(&[1, 0, 0, 0]));
        // assert_eq!(ancestors[1].state, BitVec::from_bits(&[0, 1, 0, 0]));
        // assert_eq!(ancestors[2].state, BitVec::from_bits(&[0, 0, 1, 0]));
        // assert_eq!(ancestors[3].state, BitVec::from_bits(&[0, 0, 0, 1]));
    }

    #[test]
    fn compute_multi_focal_ancestors() {
        let gen1 = BitVec::from_bits(&[0, 0, 0, 1, 1]);
        let gen2 = BitVec::from_bits(&[0, 1, 1, 0, 0]);
        let gen3 = BitVec::from_bits(&[0, 1, 1, 0, 0]);
        let gen4 = BitVec::from_bits(&[0, 0, 0, 1, 1]);

        let ag = AncestorGenerator {
            sites: vec![
                VariantSite::new(gen1, 1),
                VariantSite::new(gen2, 2),
                VariantSite::new(gen3, 3),
                VariantSite::new(gen4, 4),
            ],
        };

        let ancestors = ag.generate_ancestors();

        // TODO find out why tsinfer generates the root ancestor twice
        // TODO generate root ancestor
        assert_eq!(ancestors.len(), 3); // TODO 4

        // TODO bitvec needs Eq implementation
        // assert_eq!(ancestors[1].state, BitVec::from_bits(&[0, 1, 1, 0]));
        // assert_eq!(ancestors[2].state, BitVec::from_bits(&[1, 0, 0, 1]));
    }
}