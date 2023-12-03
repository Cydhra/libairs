use vers_vecs::{BitVec, RsVec};

/// A single variant site defined by the genotype state
#[derive(Clone, Debug)]
pub(crate) struct VariantSite {
    pub(crate) genotypes: RsVec,
    // ancestral states per sample
    pub(crate) position: usize,
    // position in the genome
    pub(crate) relative_age: f64,
}

impl VariantSite {
    pub fn new(genotypes: BitVec, position: usize) -> Self {
        let age = genotypes.count_ones() as f64 / genotypes.len() as f64;
        VariantSite {
            genotypes: RsVec::from_bit_vec(genotypes),
            position,
            relative_age: age,
        }
    }
}
