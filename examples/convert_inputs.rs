//! This binary is used to convert a VCF input file into a [`VariantData`] object and serialize it
//! into a file, so it can be used as input for the tree sequence generation.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::PathBuf;
use std::process::exit;

use clap::Parser;
use vcfire::VcfFile;

use libairs::variants::{VariantData, VariantDataBuilder};

#[derive(Parser)]
struct CliArgs {
    /// The input VCF file to convert.
    #[arg(short, long)]
    input: String,

    /// Directory where to store the output files. Optional. Defaults to the input file's directory
    #[arg(short, long)]
    output: Option<String>,

    #[arg(short, long, default_value_t = false)]
    compressed: bool,

    /// List of names of contigs to include in the output. Optional. Defaults to all contigs.
    filter: Vec<String>,
}

fn main() {
    let args = CliArgs::parse();

    let input_file = PathBuf::from(&args.input);

    let output = args.output.map_or_else(
        || PathBuf::from(input_file.parent().unwrap()),
        |name| PathBuf::from(name),
    );

    println!("reading VCF file {}", args.input);
    let variant_data = generate_variants(&args.input, args.compressed, &args.filter)
        .expect("failed to convert VCF to variant data");

    for (id, data) in variant_data {
        println!("generated {} variants for contig {}", data.len(), id);

        let output_file = output.join(format!("{}.variants", id));
        let mut output = File::create(output_file).unwrap_or_else(|error| {
            eprintln!("failed to create output file for {}: {}", id, error);
            exit(-1);
        });
        let data = serde_pickle::to_vec(&data, Default::default())
            .expect("failed to serialize variant data");
        output
            .write_all(&data)
            .expect("failed to write variant data");
    }
}

/// Read the provided VCF file and extract one variant data instance per chromosome
pub fn generate_variants(
    file: &str,
    compressed: bool,
    filter: &[String],
) -> io::Result<Vec<(String, VariantData)>> {
    let input = VcfFile::parse(file, compressed)?;

    let contigs = input
        .header
        .values
        .iter()
        .filter_map(|(name, value)| {
            if name == "contig" {
                let mut config = value[1..value.len() - 1].split(',');

                let id = config
                    .find(|s| s.starts_with("ID="))
                    .unwrap()
                    .split('=')
                    .nth(1)
                    .unwrap();
                let sequence_length = config
                    .find(|s| s.starts_with("length="))
                    .unwrap()
                    .split('=')
                    .nth(1)
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();

                Some((String::from(id), sequence_length))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    println!(
        "found {} contigs, {} of which match/es the filter.",
        contigs.len(),
        contigs
            .iter()
            .filter(|(id, _)| filter.is_empty() || filter.contains(id))
            .count()
    );

    let mut variant_data = contigs
        .iter()
        .filter_map(|(id, len)| {
            if filter.is_empty() || filter.contains(id) {
                Some((id.clone(), VariantDataBuilder::new(*len)))
            } else {
                None
            }
        })
        .collect::<HashMap<_, _>>();

    let mut counter = 0usize;

    input
        .records()?
        .map(|record| record.ok().unwrap())
        .filter(|record| filter.is_empty() || filter.contains(&record.chromosome))
        .filter(|record| {
            record.reference_bases.len() == 1
                && record.alternate_bases.len() == 1
                && record.alternate_bases[0].is_some()
                && record.alternate_bases[0].as_ref().unwrap().len() == 1
        })
        .for_each(|record| {
            let genotypes = record
                .sample_info
                .unwrap()
                .samples()
                .flat_map(|s| {
                    s.get_genotype()
                        .unwrap()
                        .split('|')
                        .map(|s| s.parse::<u8>().unwrap())
                        .collect::<Vec<_>>()
                })
                .collect();
            let position = record.position as usize;
            let reference = record.reference_bases.chars().next().unwrap();
            let alternate = record.alternate_bases[0]
                .as_ref()
                .unwrap()
                .chars()
                .next()
                .unwrap();

            // check here because of mutable reference to variant_data
            if variant_data[&record.chromosome].last_position() < Some(position) {
                variant_data
                    .get_mut(&record.chromosome)
                    .expect("filtered sample not present")
                    .add_variant_site(genotypes, position, reference, alternate);

                counter += 1;
                if counter % 10000 == 0 {
                    println!("processed {} records...", counter);
                }
            }
        });

    Ok(variant_data
        .into_iter()
        .map(|(id, v)| (id, v.finalize()))
        .collect())
}
