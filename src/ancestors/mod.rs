use std::fmt::{Debug, Display, Formatter};
use std::io;
use std::io::Write;
use std::ops::{Deref, DerefMut, Index, IndexMut};

use crate::variants::{VariantIndex, VariantSequence};
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;

mod ancestor_array;
mod generator;

pub use ancestor_array::AncestorArray;
pub use generator::AncestorGenerator;

const ANCESTRAL_STATE: u8 = 0;
const DERIVED_STATE: u8 = 1;

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
#[derive(Clone)]
pub struct AncestralSequence {
    state: VariantSequence,
    focal_sites: Vec<VariantIndex>,
    /// start of valid data in the state vector, inclusive
    start: VariantIndex,
    /// end of valid data in the state vector, exclusive
    end: VariantIndex,
    age: f64,
}

impl AncestralSequence {
    fn from_ancestral_state(len: usize, age: f64) -> Self {
        AncestralSequence {
            state: VariantSequence::from_ancestral_state(len),
            focal_sites: Vec::new(),
            start: VariantIndex(0),
            end: VariantIndex(len),
            age,
        }
    }

    /// Get the haplotype sequence for the ancestral sequence. This sequence starts at the first
    /// known site and ends at the last known site. Where this sequence is located in the genome
    /// is defined by [`start`] and [`end`] respectively.
    /// This means that the indices of the haplotype sequence do not correspond to the indices of
    /// the genome.
    ///
    /// [`start`]: Self::start
    /// [`end`]: Self::end
    pub fn haplotype(&self) -> &[u8] {
        &self.state[self.start..self.end]
    }

    /// Get an enumerated iterator over the haplotype sequence. This sequence starts at the first
    /// known site and ends at the last known site. Each returned element also contains the index
    /// of the site in the genome (regarding the genome's variant site vector, the actual genome
    /// position will differ from this).
    pub(crate) fn site_iter(
        &self,
    ) -> impl Iterator<Item = (VariantIndex, &'_ u8)> + DoubleEndedIterator + '_ {
        self.state
            .iter()
            .enumerate()
            .skip(self.start.0)
            .take(self.start.get_variant_distance(self.end))
            .map(|(idx, b)| (VariantIndex(idx), b))
    }

    /// Get the position of the first known site in the genome (inclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub(crate) fn start(&self) -> VariantIndex {
        self.start
    }

    /// Get the position of the last known site in the genome (exclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub(crate) fn end(&self) -> VariantIndex {
        self.end
    }

    /// Get the inferred relative age of the ancestral sequence. This is derived from the inferred
    /// relative age of the focal sites that were used to infer this ancestor.
    pub fn relative_age(&self) -> f64 {
        self.age
    }

    /// Get the length of the ancestral sequence. Only known sites are considered, so the length
    /// might be shorter than the length of the underlying DNA sequence.
    pub fn len(&self) -> usize {
        self.end.get_variant_distance(self.start)
    }

    /// Dump the ancestral sequence into a text file for the testing suite.
    pub fn export(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{start}\t{end}\t{age}\t{focal_sites:?}\t",
            start = self.start.0,
            end = self.end.0,
            age = self.age,
            focal_sites = self.focal_sites.iter().map(|s| s.0).collect::<Vec<_>>(),
        ))?;

        for b in &self.state[self.start..self.end] {
            writer.write_fmt(format_args!("{}", b))?;
        }

        writer.write_fmt(format_args!("\n"))?;
        Ok(())
    }
}

impl Debug for AncestralSequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("AncestralSequence { ")?;
        f.write_fmt(format_args!("focal_sites={:?},\t", self.focal_sites))?;

        f.write_str("genotype=[ ")?;
        let mut iter = self.state.iter();
        let mut idx = 0;
        while let Some(b) = iter.next() {
            if idx < self.start.0 || idx >= self.end.0 {
                f.write_str("-")?;
            } else {
                f.write_fmt(format_args!("{}", b))?;
            }

            if iter.len() > 0 {
                f.write_str(", ")?;
            }

            idx += 1;
        }
        f.write_str(" ] }")?;
        Ok(())
    }
}

// TODO can we get rid of this implementation to force the new-type?
impl Index<usize> for AncestralSequence {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.state.deref()[index]
    }
}

impl IndexMut<usize> for AncestralSequence {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.state.deref_mut()[index]
    }
}

impl Index<VariantIndex> for AncestralSequence {
    type Output = u8;

    fn index(&self, index: VariantIndex) -> &Self::Output {
        &self.state[index]
    }
}

impl IndexMut<VariantIndex> for AncestralSequence {
    fn index_mut(&mut self, index: VariantIndex) -> &mut Self::Output {
        &mut self.state[index]
    }
}

/// An index into the ancestor array which uniquely identifies an ancestor
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) struct Ancestor(pub(crate) usize);

impl Display for Ancestor {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}
