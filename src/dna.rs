use std::fmt::{Display, Formatter};

/// A single variant site defined by the genotype state.
#[derive(Clone, Debug)]
pub struct VariantSite {
    // todo hide most fields
    pub(crate) genotypes: Vec<u8>,
    // ancestral states per sample
    pub position: SequencePosition,
    // position in the genome
    pub(crate) relative_age: f64,
    // whether the site is bi-allelic or not
    pub(crate) is_biallelic: bool,
    // whether the site is a singleton or not
    pub(crate) is_singleton: bool,
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
    /// - `position` the position of the site in the genome.
    pub fn new(genotypes: Vec<u8>, position: usize) -> Self {
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
            position: SequencePosition(position),
            relative_age: age,
            is_biallelic,
            is_singleton,
        }
    }
}

/// A position in a DNA sequence. This newtype ensures that sequence positions and variant indices (indices into the
/// variant site vector) aren't mixed up.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SequencePosition(usize); // TODO do not make this public

impl SequencePosition {
    /// Create a new sequence position from a usize.
    pub fn from_usize(position: usize) -> Self {
        Self(position)
    }

    #[inline]
    pub fn from_vec(positions: Vec<usize>) -> Vec<Self> {
        positions.into_iter().map(|p| Self::from_usize(p)).collect()
    }

    /// Get the underlying usize value of the sequence position.
    pub fn unwrap(&self) -> usize {
        self.0
    }
}

impl Display for SequencePosition {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
