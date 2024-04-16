use crate::variants::sequence::MutationState;
use crate::variants::{SequencePosition, VariantData, VariantSite};

/// A builder for [`VariantData`] instances
#[derive(Clone, Debug)]
pub struct VariantDataBuilder {
    sites: Vec<VariantSite>,
    positions: Vec<SequencePosition>,
    num_samples: usize,
    last_position: Option<SequencePosition>,
    sequence_length: SequencePosition,
}

impl VariantDataBuilder {
    /// A new empty builder for a genome with given sequence length
    pub fn new(sequence_length: usize) -> Self {
        Self {
            sites: Vec::default(),
            positions: Vec::default(),
            num_samples: 0,
            last_position: None,
            sequence_length: SequencePosition::from_usize(sequence_length),
        }
    }

    /// Construct a new builder from an iterator over variant sites and their sequence positions.
    ///
    /// # Parameters
    /// - `sequence_length` the length of the reference genome
    /// - `iter` an iterator over triples of variant sites and their sequence positions in
    /// increasing order of sequence position, as well as the derived state in FASTA notation
    pub fn from_iter<I>(sequence_length: usize, iter: I) -> Self
    where
        I: IntoIterator<Item = (Vec<MutationState>, usize, char, char)>,
    {
        let mut builder = Self::new(sequence_length);
        for (site, position, ancestral_state, derived_state) in iter {
            builder.add_variant_site(site, position, ancestral_state, derived_state);
        }
        builder
    }

    /// Add a variant site to the variant data. If the variant site cannot be used for inference
    /// (for example, if it is a singleton mutation), it will be ignored.
    /// This behavior will likely change in the future, so the algorithm can at least encode the
    /// additional information into the final tree sequence, even if not used for inference.
    /// Mutation sites must be added in increasing order of sequence position.
    ///
    /// # Parameters
    /// - `state` a vector of [`MutationState`] values that indicate the variant call for each
    /// sample. The number of samples is inferred from the first state vector added to the builder.
    /// The state vector must not be empty.
    /// - `sequence_position` the position of the mutation site in the reference genome
    /// - `derived_state` the derived state of the mutation site in FASTA notation
    ///
    /// # Panics
    /// - if the state vector is empty
    /// - if the state vector has a different number of entries than previously added state vectors
    /// - if a position is added that is not greater than the last position added
    pub fn add_variant_site(
        &mut self,
        state: Vec<MutationState>,
        sequence_position: usize,
        ancestral_state: char,
        derived_state: char,
    ) {
        // TODO store invalid variant sites somewhere for encoding in final tree sequence
        assert!(!state.is_empty());
        assert!(self.num_samples == 0 || state.len() == self.num_samples);

        if self.num_samples == 0 {
            self.num_samples = state.len();
        }

        let sequence_position = SequencePosition::from_usize(sequence_position);

        if self.last_position < Some(sequence_position) {
            self.last_position = Some(sequence_position);
        } else {
            panic!("Variant sites must be added in increasing order of sequence position");
        }

        let variant_site =
            VariantSite::new(state, sequence_position, ancestral_state, derived_state);
        if Self::is_valid_site(&variant_site) {
            self.positions.push(sequence_position);
            self.sites.push(variant_site)
        }
    }

    /// Generate [`VariantData`] from the current builder state.
    pub fn finalize(self) -> VariantData {
        VariantData::new(
            self.sites,
            self.positions,
            self.sequence_length,
            self.num_samples,
        )
    }

    /// Whether a site is valid for the algorithm. A site is valid if it is not a singleton and
    /// is biallelic.
    fn is_valid_site(site: &VariantSite) -> bool {
        !site.is_singleton && site.is_biallelic
    }
}
