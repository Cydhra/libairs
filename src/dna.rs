use std::fmt::{Debug, Formatter, Write};
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
        self.sequence.get_bits_unchecked(index * self.encoding as usize, self.encoding as usize)
    }

    pub(crate) fn get(&self, index: usize) -> Option<u64> {
        self.sequence.get_bits(index * self.encoding as usize, self.encoding as usize)
    }
}

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
// TODO: is this enough? If it isn't, this must be replaced with a packed sequence coding for whatever
//  state we need, ideally with a optimization because we expect most of the entries to reference the
//  ancestral state
pub(crate) struct AncestralSequence {
    pub(crate) state: BitVec,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl AncestralSequence {
    pub(crate) fn from_ancestral_state(len: usize) -> Self {
        AncestralSequence {
            state: BitVec::from_zeros(len),
            start: 0,
            end: 0,
        }
    }

    pub(crate) fn set_unchecked(&mut self, index: usize, value: u64) {
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