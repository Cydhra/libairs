use crate::ancestors::AncestorGenerator;
use crate::dna::VariantSite;
use std::io;
use vcfire::VcfFile;

pub fn from_vcf(file: &str, compressed: bool) -> io::Result<AncestorGenerator> {
    let input = VcfFile::parse(file, compressed)?;
    Ok(AncestorGenerator::from_iter(input.records()?.map(
        |record| {
            VariantSite::new(
                record?
                    .sample_info
                    .iter()
                    .flat_map(|sample_info| {
                        sample_info
                            .samples()
                            .map(|s| s.get_genotype().unwrap().split('|'))
                    })
                    .collect(),
                record?.position as usize,
            )
        },
    )))
}
