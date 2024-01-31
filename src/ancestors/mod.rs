use std::fmt::{Debug, Formatter};
use std::io;
use std::io::Write;
use std::ops::{Index, IndexMut};

use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;

mod ancestor_array;
mod generator;

pub(crate) use ancestor_array::AncestorArray;
pub use generator::AncestorGenerator;

const ANCESTRAL_STATE: u8 = 0;
const DERIVED_STATE: u8 = 1;

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
#[derive(Clone)]
pub struct AncestralSequence {
    state: Vec<u8>,
    focal_sites: Vec<usize>,
    /// start of valid data in the state vector, inclusive
    start: usize,
    /// end of valid data in the state vector, exclusive
    end: usize,
    age: f64,
}

impl AncestralSequence {
    // TODO shouldnt be public. instead a builder should be used, so the meta data can be calculated
    pub fn from_ancestral_state(len: usize, age: f64) -> Self {
        AncestralSequence {
            state: vec![0u8; len],
            focal_sites: Vec::new(),
            start: 0,
            end: 0,
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
    pub fn site_iter(&self) -> impl Iterator<Item = (usize, &'_ u8)> + DoubleEndedIterator + '_ {
        self.state
            .iter()
            .enumerate()
            .skip(self.start)
            .take(self.end - self.start)
    }

    /// Get the position of the first known site in the genome (inclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub fn start(&self) -> usize {
        self.start
    }

    /// Get the position of the last known site in the genome (exclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub fn end(&self) -> usize {
        self.end
    }

    /// Get the set of focal sites that were used to generate this ancestral sequence.
    pub fn focal_sites(&self) -> &[usize] {
        &self.focal_sites
    }

    /// Get the inferred relative age of the ancestral sequence. This is derived from the inferred
    /// relative age of the focal sites that were used to infer this ancestor.
    pub fn relative_age(&self) -> f64 {
        self.age
    }

    /// Get the length of the ancestral sequence. Only known sites are considered, so the length
    /// might be shorter than the length of the underlying DNA sequence.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Dump the ancestral sequence into a text file for the testing suite.
    pub fn export(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{start}\t{end}\t{age}\t{focal_sites:?}\t",
            start = self.start,
            end = self.end,
            age = self.age,
            focal_sites = self.focal_sites
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
            if idx < self.start || idx >= self.end {
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

impl Index<usize> for AncestralSequence {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.state[index]
    }
}

impl IndexMut<usize> for AncestralSequence {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.state[index]
    }
}

/// An index into the [`AncestorArray`]
///
/// [`AncestorArray`]: AncestorArray
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) struct VariantIndex(usize);

impl VariantIndex {
    #[cfg(test)]
    pub fn from_usize(index: usize) -> Self {
        Self(index)
    }

    /// Get the next variant index after this one
    pub(crate) fn next(&self) -> Self {
        Self(self.0 + 1)
    }
}

/// An index into the ancestor array which uniquely identifies an ancestor
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) struct Ancestor(pub(crate) usize);
