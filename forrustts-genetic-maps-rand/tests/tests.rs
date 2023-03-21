use forrustts_genetic_maps::GeneticMapBuilder;
use forrustts_genetic_maps_rand::GeneticMap;

#[test]
fn test_build_empty() {
    let builder = GeneticMapBuilder::default();
    assert!(GeneticMap::new_from_builder(builder).is_some());
}

#[test]
fn test_build_non_empty() {
    unimplemented!("test not written");
}

#[test]
fn test_generate_breakpoints() {
    unimplemented!("test not written");
}
