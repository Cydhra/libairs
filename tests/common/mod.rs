use libairs::ancestors::AncestorGenerator;
use libairs::variants::VariantDataBuilder;

pub fn create_ancestor_generator<const S: usize>(
    sequence_len: usize,
    variant_sites: &[[u8; S]],
) -> AncestorGenerator {
    let variant_data = VariantDataBuilder::from_iter(
        sequence_len,
        variant_sites
            .iter()
            .enumerate()
            .map(|(i, site)| (site.to_vec(), i + 1)),
    )
    .finalize();

    AncestorGenerator::from_variant_data(variant_data)
}
