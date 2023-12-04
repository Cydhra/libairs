use vers_vecs::BitVec;

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