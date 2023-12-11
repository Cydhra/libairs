use crate::ancestors::AncestralSequence;

struct TreeSequenceGenerator {
    recombination_probabilities: Vec<f64>,
    mismatch_probabilities: Vec<f64>,
}

impl TreeSequenceGenerator {
    pub fn new(recombination_rate: f64, mismatch_rate: f64, variant_positions: Vec<usize>) -> Self {
        // TODO those values are only tsinfer's defaults, the actual calculation of the values
        //  works differently

        let mut recombination_probabilities = vec![recombination_rate; variant_positions.len()];
        recombination_probabilities[0] = 0f64;

        let mismatch_probabilities = vec![mismatch_rate; variant_positions.len()];

        Self {
            recombination_probabilities,
            mismatch_probabilities,
        }
    }

    /// For a given [`AncestralSequence`] and a set of ancestor sequences, calculate the most likely
    /// copying path within the LS model using the Viterbi algorithm.
    // TODO the ancestors are probably required to be packed sequences instead of ancestral sequences
    // TODO ancestors are required to be TSNodeSequences, but those dont work yet
    fn find_hidden_path(
        &self,
        candidate: &AncestralSequence,
        ancestors: &[AncestralSequence],
    ) -> Vec<(usize, usize, usize)> {
        let mut likelihoods = vec![1f64; ancestors.len()];
        let mut recombination_points = vec![vec![false; ancestors.len()]; candidate.len()];
        let mut max_likelihoods = vec![0; candidate.len()];
        let candidate_start = candidate.start();

        // println!("=================================================================");
        // println!("CALCULATING FOR ANCESTOR {:?}", candidate);
        for (site, &state) in candidate.site_iter() {
            let mut max_site_likelihood = -1f64;
            let mut max_site_likelihood_ancestor: Option<usize> = None;

            let rho = self.recombination_probabilities[site];
            let mu = self.mismatch_probabilities[site];
            let k = ancestors.len() as f64;
            let num_alleles = 2f64; // TODO we might not want to hard-code this

            // todo only process sites that are actually valid in the ancestor
            for (ancestor, ancestral_sequence) in ancestors.iter().enumerate() {
                // println!("p_last: {}", likelihoods[ancestor]);
                let prob_no_recomb = likelihoods[ancestor] * (1f64 - rho - rho / k);
                // println!("prob_no_recomb: {}", prob_no_recomb);
                let prob_recomb = rho / k;
                // println!("prob_recomb: {}", prob_recomb);

                let pt = if prob_no_recomb > prob_recomb {
                    prob_no_recomb
                } else {
                    recombination_points[site - candidate_start][ancestor] = true;
                    prob_recomb
                };

                // println!("mu: {}", mu);
                // TODO ancestral sequences return garbage data on out-of-bounds queries, catch that!
                let pe = if state == ancestral_sequence[site] {
                    1f64 - (num_alleles - 1f64) * mu
                } else {
                    mu
                };
                // println!("pe: {}", pe);

                likelihoods[ancestor] = pt * pe;
                // println!(
                //     "updated likelihood at {} to {}",
                //     ancestor, likelihoods[ancestor]
                // );
                if likelihoods[ancestor] > max_site_likelihood {
                    max_site_likelihood = likelihoods[ancestor];
                    max_site_likelihood_ancestor = Some(ancestor);
                }
            }

            // Apparently a measure to maintain numerical stability
            for (ancestor, _) in ancestors.iter().enumerate() {
                likelihoods[ancestor] /= max_site_likelihood;
            }

            // println!(
            //     "site {site} maximum likelihood: {max_site_likelihood} selecting ancestor {}",
            //     max_site_likelihood_ancestor.unwrap()
            // );
            max_likelihoods[site - candidate_start] =
                max_site_likelihood_ancestor.expect("no max likelihood calculated");
        }

        let mut nodes = Vec::new();
        let mut ancestor_index = max_likelihoods[candidate.end() - 1 - candidate_start];
        let mut ancestor_coverage_end = candidate.end();

        for (site, &state) in candidate.site_iter().rev().skip(1) {
            if recombination_points[site - candidate_start][ancestor_index] {
                nodes.push((ancestor_index, site, ancestor_coverage_end));
                assert_ne!(ancestor_index, max_likelihoods[site - 1 - candidate_start]);
                ancestor_index = max_likelihoods[site - 1 - candidate_start];
                ancestor_coverage_end = site;
            }
        }
        nodes.push((ancestor_index, 0, ancestor_coverage_end));

        nodes.reverse();
        nodes
    }

    pub fn generate_tree_sequence(
        &self,
        ancestors: Vec<AncestralSequence>,
        ancestral_state: AncestralSequence,
    ) {
        debug_assert!(ancestors
            .windows(2)
            .all(|window| window[0].relative_age() >= window[1].relative_age()));

        let mut tableau = Vec::with_capacity(ancestors.len());
        tableau.push(ancestral_state.clone());

        let mut current_age = f64::INFINITY;
        let mut current_age_set = Vec::new();
        let mut tree: Vec<(Vec<u8>, Vec<(usize, usize, usize)>)> = Vec::new();
        for ancestor in ancestors {
            if ancestor.relative_age() < current_age {
                tableau.append(&mut current_age_set);
                current_age = ancestor.relative_age();
            }

            tree.push((
                Vec::from(ancestor.haplotype().clone()),
                self.find_hidden_path(&ancestor, &tableau),
            ));
            current_age_set.push(ancestor);
        }

        // println!("=================================================================");
        // println!("FINAL TREE");
        // println!(
        //     "ancestor: {:?}: \t{:?}",
        //     ancestral_state.haplotype(),
        //     vec![(-1, 0, 4)]
        // );
        // for (ancestor, path) in tree {
        //     println!("ancestor: {:?}: \t{:?}", ancestor, path);
        // }
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::{AncestorGenerator, AncestralSequence};
    use crate::dna::VariantSite;
    use crate::ts::TreeSequenceGenerator;

    #[test]
    fn trivial_tree_test() {
        let site1 = vec![0, 0, 0, 1, 1, 1];
        let site2 = vec![0, 1, 1, 0, 0, 0];
        let site3 = vec![0, 1, 1, 0, 0, 0];
        let site4 = vec![0, 0, 0, 1, 1, 1];
        let site5 = vec![0, 1, 0, 0, 0, 1];

        let ag = AncestorGenerator::from_iter(
            vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site3, 3),
                VariantSite::new(site4, 4),
                VariantSite::new(site5, 5),
            ]
                .into_iter(),
        );

        let mut ancestors = ag.generate_ancestors();

        let ancestor_matcher = TreeSequenceGenerator::new(1e-2, 1e-20, vec![1, 2, 3, 4, 5, 6]);

        ancestors.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());
        for ancestor in &ancestors {
            println!(
                "ancestral sequence: {:?} has age {}",
                ancestor,
                ancestor.relative_age()
            );
        }

        let mut ancestral_state = AncestralSequence::from_ancestral_state(5, 2.0f64);
        ancestral_state.end = 5;
        ancestor_matcher.generate_tree_sequence(ancestors, ancestral_state);

        // TODO verify results
    }

    #[test]
    fn trivial_recombinant_test() {
        let site1 = vec![0, 0, 0, 1, 1, 1];
        let site2 = vec![0, 0, 0, 1, 1, 1];
        let site4 = vec![0, 1, 0, 0, 0, 1];
        let site5 = vec![0, 0, 0, 1, 1, 0];
        let site6 = vec![0, 1, 1, 0, 0, 1];
        let site7 = vec![0, 1, 1, 0, 0, 1];

        let ag = AncestorGenerator::from_iter(
            vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site4, 4),
                VariantSite::new(site5, 5),
                VariantSite::new(site6, 6),
                VariantSite::new(site7, 7),
            ]
                .into_iter(),
        );

        let mut ancestors = ag.generate_ancestors();
        let ancestor_matcher = TreeSequenceGenerator::new(1e-2, 1e-20, vec![1, 2, 4, 5, 6, 7]);

        ancestors.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());
        for ancestor in &ancestors {
            println!(
                "ancestral sequence: {:?} has age {}",
                ancestor,
                ancestor.relative_age()
            );
        }

        let mut ancestral_state = AncestralSequence::from_ancestral_state(6, 2.0f64);
        ancestral_state.end = 6;
        ancestor_matcher.generate_tree_sequence(ancestors, ancestral_state);

        // TODO verify results
    }
}
