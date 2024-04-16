use crate::ancestors::AncestorGenerator;
use crate::variants::VariantDataBuilder;
use std::io;
use vcfire::VcfFile;

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
