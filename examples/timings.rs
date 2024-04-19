//! This example is used by a python module for benchmarking airs. It offers a CLI to control its
//! behavior and measures time itself to make it easier to exclude parsing times for data.

use std::io;
use std::path::{Path, PathBuf};
use std::process::exit;
use std::time::{Duration, Instant};

use clap::{Args, Parser, Subcommand, ValueEnum};

use libairs::ancestors::AncestorGenerator;
use libairs::variants::{VariantData, VariantDataBuilder};

#[derive(Parser)]
#[command(version, arg_required_else_help = true)]
struct CliArgs {
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    num_threads: u16,

    #[command(subcommand)]
    command: Action,
}

#[derive(Subcommand)]
enum Action {
    GenerateAncestors {
        #[command(flatten)]
        data_source: Input,

        #[arg(short, long)]
        sequence_length: usize,

        #[arg(short, long)]
        output: Option<String>,
    },
    MatchAncestors {},
    MatchSamples {},
    Infer {
        #[command(flatten)]
        data_source: Input,

        #[arg(short, long)]
        sequence_length: usize,
    },
}

#[derive(Args)]
#[group(required = true)]
struct Input {
    #[arg(long = "type", value_enum)]
    input_type: InputType,

    #[arg(short = 'i', long)]
    path: String,

    #[arg(long, default_value_t = false)]
    compressed: bool,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum InputType {
    /// Import from a VCF file (proprietary vcf subset)
    Vcf,

    /// Import from a file containing sample data directly printed by a python script
    Python,
}

fn main() {
    let args = CliArgs::parse();

    if args.num_threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads as usize)
            .build_global()
            .unwrap();
    }

    match args.command {
        Action::GenerateAncestors {
            data_source,
            sequence_length,
            output,
        } => {
            let input_path = data_source.path.clone();

            let ag = parse_input(data_source, sequence_length).unwrap_or_else(|error| {
                eprintln!("could not parse input data: {}", error);
                exit(-1);
            });

            let start = Instant::now();
            let ancestors = ag.generate_ancestors();
            let end = start.elapsed();
            println!("generated {} ancestors in {:?}", ancestors.len(), end);

            let output = if let Some(path) = output {
                path
            } else {
                String::from(
                    PathBuf::from(&input_path).parent().unwrap_or_else(|| {
                        eprintln!("cannot export ancestors because input file has no parent and no output path was provided.");
                        exit(-1);
                }).to_str().unwrap())
            };

            ancestors
                .export_ancestors(Path::new(&output))
                .unwrap_or_else(|error| {
                    eprintln!("failed to export ancestors to {}: {}", output, error)
                })
        }
        Action::MatchAncestors { .. } => {}
        Action::MatchSamples { .. } => {}
        Action::Infer {
            data_source,
            sequence_length,
        } => {
            let ag = parse_input(data_source, sequence_length).unwrap_or_else(|error| {
                eprintln!("could not parse input data: {}", error);
                exit(-1);
            });

            let mut total = Duration::new(0, 0);

            let start = Instant::now();
            let ancestors = ag.generate_ancestors();
            let end = start.elapsed();
            total += end;
            println!("generated {} ancestors in {:?}", ancestors.len(), end);

            let start = Instant::now();
            let mut ancestor_matcher =
                libairs::ts::ViterbiMatcher::new(ancestors, 1e-2, 1e-20, true, 40);
            ancestor_matcher.match_ancestors();
            let end = start.elapsed();
            total += end;
            println!("matched ancestors in {:?}", end);

            let start = Instant::now();
            ancestor_matcher.match_samples();
            let end = start.elapsed();
            total += end;
            println!("matched samples in {:?}", end);

            println!("total time: {:?}", total);
        }
    }
}

fn parse_input(data_source: Input, sequence_length: usize) -> io::Result<AncestorGenerator> {
    match data_source.input_type {
        InputType::Vcf => libairs::convenience::from_vcf(
            &data_source.path,
            data_source.compressed,
            sequence_length,
        ),
        InputType::Python => import_python_data(&data_source.path, sequence_length),
    }
}

/// Import a custom file with the following layout:
/// ```
/// num_samples
/// 0 1 0 1 repeat #num_samples... 0 position ancestral_state derived_state
/// repeat #num_variants often...
/// ```
fn import_python_data(path: &String, sequence_length: usize) -> io::Result<AncestorGenerator> {
    let data = std::fs::read_to_string(path)?;
    let mut lines = data.lines();
    let num_samples = lines.next()?.parse::<usize>().unwrap();
    let mut builder = VariantDataBuilder::new(sequence_length);
    for line in lines {
        let mut parts = line.split_whitespace();
        let mut states = Vec::with_capacity(num_samples);
        for _ in 0..num_samples {
            states.push(parts.next()?.parse::<u8>().unwrap());
        }
        let position = parts.next()?.parse::<usize>().unwrap();
        let ancestral_state = parts.next()?.chars().next().unwrap();
        let derived_state = parts.next()?.chars().next().unwrap();
        builder.add_variant_site(states, position, ancestral_state, derived_state);
    }
    Ok(AncestorGenerator::new(builder.finalize()))
}
