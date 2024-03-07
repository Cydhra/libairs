use crate::ancestors::AncestorGenerator;
use crate::variants::VariantSite;
use std::io;
use vcfire::VcfFile;

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
            record.position as usize,
        )
    }));
    Ok(generator)
}
