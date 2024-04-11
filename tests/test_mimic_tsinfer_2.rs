//! tests whether airs mimics tsinfer's output

mod common;

#[test]
fn test_mimic_tsinfer_2() {
    let sites = [
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ];

    let ag = common::create_ancestor_generator(6, &sites);
    let ancestors = common::generate_ancestors(ag);
    let ts = common::match_ancestors(ancestors);

    // test whether the tsinfer output is generated

    assert_eq!(ts[1].edges().len(), 1);
    assert_eq!(ts[1].edges()[0].parent, 0);

    assert_eq!(ts[2].edges().len(), 1);
    assert_eq!(ts[2].edges()[0].parent, 1);

    assert_eq!(ts[3].edges().len(), 1);
    assert_eq!(ts[3].edges()[0].parent, 2);

    assert_eq!(ts[4].edges().len(), 1);
    assert_eq!(ts[4].edges()[0].parent, 3);

    assert_eq!(ts[5].edges().len(), 1);
    assert_eq!(ts[5].edges()[0].parent, 3);

    assert_eq!(ts[6].edges().len(), 3);
    assert!(ts[6]
        .edges()
        .iter()
        .any(|ni| ni.parent == 5 || ni.parent == 4));
    assert!(ts[6].edges().iter().any(|ni| ni.parent == 2));
    assert!(ts[6].edges().iter().any(|ni| ni.parent == 0));
}
