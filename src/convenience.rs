use crate::ancestors::AncestorGenerator;
use crate::variants::{SequencePosition, VariantSite};
use std::io;
use vcfire::VcfFile;

// TODO change this implementation to use the builder interface and VariantData instead
pub fn from_vcf(file: &str, compressed: bool) -> io::Result<AncestorGenerator> {
    let input = VcfFile::parse(file, compressed)?;
    let generator = AncestorGenerator::from_iter(input.records()?.map(|record| {
        let record = record.ok().unwrap();
        VariantSite::new(
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
            SequencePosition::from_usize(record.position as usize),
        )
    }));
    Ok(generator)
}
