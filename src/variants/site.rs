use crate::variants::SequencePosition;

/// A single variant site defined by the genotype state.
#[derive(Clone, Debug)]
pub struct VariantSite {
    pub(crate) genotypes: Vec<u8>,
    // position in the genome
    pub(crate) position: SequencePosition,
    // putative age of this site
    pub(crate) relative_age: f64,
    // whether the site is bi-allelic
    pub(crate) is_biallelic: bool,
    // whether the site is a singleton
    pub(crate) is_singleton: bool,
    // the genomic states in FASTA notation
    pub(crate) ancestral_state: char,
    pub(crate) derived_state: char,
}

impl VariantSite {
    /// Create a new variant site from a vector of genotypes and its position in the genome. At the
    /// moment, only biallelic sites are supported, meaning any entry in the vector must be either
    /// 0 or 1. Entries that are greater than 1 will be treated as 1.
    ///
    /// # Parameters
    /// - `genotypes` a vector if `n` integers, each defining which allele is present at the site.
    /// 0 is the reference allele, 1 is the derived allele. Any value greater than 1 will be
    /// treated as 1.
    /// - `position` the [`SequencePosition`] of the site in the genome.
    pub fn new(
        genotypes: Vec<u8>,
        position: SequencePosition,
        ancestral_state: char,
        derived_state: char,
    ) -> Self {
        let mut derived_sites = 0;
        let mut highest_state = 0;
        genotypes.iter().for_each(|&state| {
            if state > 0 {
                derived_sites += 1;
            }
            if state > highest_state {
                highest_state = state;
            }
        });

        let age = derived_sites as f64 / genotypes.len() as f64;
        // TODO technically this is not correct, since there could be a derived allele that isn't
        //  represented with a 1.
        let is_biallelic = highest_state == 1;
        let is_singleton = derived_sites == 1;
        VariantSite {
            genotypes,
            position,
            relative_age: age,
            is_biallelic,
            is_singleton,
            ancestral_state,
            derived_state,
        }
    }
}
