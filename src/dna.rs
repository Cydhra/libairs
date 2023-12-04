/// A single variant site defined by the genotype state
#[derive(Clone, Debug)]
pub(crate) struct VariantSite {
    pub(crate) genotypes: Vec<u8>,
    // ancestral states per sample
    pub(crate) position: usize,
    // position in the genome
    pub(crate) relative_age: f64,
}

impl VariantSite {
    pub fn new(genotypes: Vec<u8>, position: usize) -> Self {
        let age = genotypes.iter().filter(|&s| *s > 0).count() as f64 / genotypes.len() as f64;
        VariantSite {
            genotypes,
            position,
            relative_age: age,
        }
    }
}
