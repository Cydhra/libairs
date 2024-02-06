//! This example shows how to generate a tree sequence from ancestors from a VCF file.
//! The actual sample data is not included in the tree sequence.
//! The tree sequence is exported to a file matching the VCF file but with the extension `.trees`.
//! It can be imported into tskit using [`tskit.load_text()`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.load_text).

use libairs::dna::SequencePosition;
use std::env;
use std::path::PathBuf;
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
    let sequence_length = SequencePosition::from_usize(
        contig_config[1..contig_config.len() - 1]
            .split(',')
            .find(|s| s.starts_with("length="))
            .unwrap()
            .split('=')
            .nth(1)
            .unwrap()
            .parse::<usize>()
            .unwrap(),
    );

    let mut target_file = PathBuf::from(&vcf);

    let ancestor_generator = libairs::convenience::from_vcf(&vcf, compressed).unwrap();
    let ancestors = ancestor_generator.generate_ancestors(sequence_length);
    let matcher = libairs::ts::TreeSequenceGenerator::new(
        ancestors,
        sequence_length,
        1e-2,
        1e-20,
        ancestor_generator.variant_positions(),
    );
    target_file.pop();
    matcher
        .export_ancestors(&target_file)
        .expect("failed to export ancestors");

    let tree_sequence = matcher.generate_tree_sequence();

    tree_sequence
        .tskit_export(&target_file)
        .expect("failed to export tree sequence");

    ancestor_generator
        .tskit_export_sites(&target_file.as_path())
        .expect("failed to export sites");
}
