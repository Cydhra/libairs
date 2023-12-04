use crate::ancestors::AncestralSequence;

struct TreeSequenceGenerator {}

impl TreeSequenceGenerator {
    /// For a given [`AncestralSequence`] and a set of ancestor sequences, calculate the most likely
    /// copying path within the LS model using the Viterbi algorithm.
    // TODO the ancestors are probably required to be packed sequences instead of ancestral sequences
    fn find_hidden_path(&self, candidate: AncestralSequence, ancestors: &[AncestralSequence]) {}

    pub fn generate_tree_sequence(&self) {}
}