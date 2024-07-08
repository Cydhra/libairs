//! This example shows how to generate a tree sequence from ancestors from a VCF file.
//! The actual sample data is not included in the tree sequence.
//! The tree sequence is exported to a file matching the VCF file but with the extension `.trees`.
//! It can be imported into tskit using [`tskit.load_text()`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.load_text).

use libairs::ancestors::AncestorGenerator;
use libairs::ts::ViterbiMatcher;
use libairs::variants::VariantDataBuilder;
use std::path::PathBuf;
use std::{env, io};
use vcfire::VcfFile;

fn main() {
    if env::args().len() < 2 {
        println!(
            "usage: {} <vcf file> [<compressed>]",
            env::args().nth(0).unwrap()
        );
        return;
    }
    let vcf = env::args().nth(1).unwrap();
    let compressed = if let Some(c) = env::args().nth(2) {
        if let Ok(c) = c.parse::<bool>() {
            c
        } else {
            println!("invalid value for <compressed> (must be true or false)");
            return;
        }
    } else {
        false
    };

    let vcf_header = VcfFile::parse(&vcf, compressed).unwrap().header;
    let contig_config = vcf_header
        .values
        .iter()
        .find(|(key, _)| key == "contig")
        .unwrap()
        .1
        .clone();
    let sequence_length = contig_config[1..contig_config.len() - 1]
        .split(',')
        .find(|s| s.starts_with("length="))
        .unwrap()
        .split('=')
        .nth(1)
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let mut target_file = PathBuf::from(&vcf);

    let ancestor_generator = from_vcf(&vcf, compressed, sequence_length).unwrap();

    target_file.pop();
    ancestor_generator
        .tskit_export_sites(&target_file.as_path())
        .expect("failed to export sites");

    let ancestors = ancestor_generator.generate_ancestors();
    ancestors
        .export_ancestors(&target_file)
        .expect("failed to export ancestors");

    let mut ancestor_matcher = ViterbiMatcher::new(ancestors, 1e-2, 1e-20, false, 1);
    ancestor_matcher.match_ancestors();
    ancestor_matcher.match_samples();

    let tree_sequence = ancestor_matcher.get_tree_sequence();

    tree_sequence
        .tskit_export(&target_file)
        .expect("failed to export tree sequence");
}

/// Convert a synthetically generated VCF file from msprime that contains a single contig into a
/// [`AncestorGenerator`] instance
pub fn from_vcf(
    file: &str,
    compressed: bool,
    sequence_length: usize,
) -> io::Result<AncestorGenerator> {
    let input = VcfFile::parse(file, compressed)?;
    let variant_data = VariantDataBuilder::from_iter(
        sequence_length,
        input
            .records()?
            .map(|record| record.ok().unwrap())
            .filter(|record| record.reference_bases.len() == 1 && record.alternate_bases.len() == 1)
            .map(|record| {
                (
                    record
                        .sample_info
                        .iter()
                        .flat_map(|sample_info| {
                            sample_info.samples().flat_map(|s| {
                                s.get_genotype()
                                    .unwrap()
                                    .split('|')
                                    .map(|s| s.parse::<u8>().unwrap())
                                    .collect::<Vec<_>>()
                            })
                        })
                        .collect(),
                    record.position as usize,
                    record.reference_bases.chars().next().unwrap(),
                    record.alternate_bases[0]
                        .as_ref()
                        .unwrap()
                        .chars()
                        .next()
                        .unwrap(),
                )
            }),
    )
    .finalize();

    Ok(AncestorGenerator::from_variant_data(variant_data))
}
