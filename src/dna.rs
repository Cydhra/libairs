use vers_vecs::{BitVec, RsVec};

// TODO benchmark if bio-seq is fast enough to replace this with a fully engineered version

/// A packed nucleotide sequence. Supports any number of bits for packing, meaning this type can pack
/// raw read data in 2 bits per base, and data containing ambiguity codes in larger strides.
pub(crate) struct PackedSequence {
    sequence: BitVec,
    // packed sequence
    encoding: u8, // how many bits are required per nucleotide
}

impl PackedSequence {
    pub(crate) fn get_unchecked(&self, index: usize) -> u64 {
        self.sequence
            .get_bits_unchecked(index * self.encoding as usize, self.encoding as usize)
    }

    pub(crate) fn get(&self, index: usize) -> Option<u64> {
        self.sequence
            .get_bits(index * self.encoding as usize, self.encoding as usize)
    }
}

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
