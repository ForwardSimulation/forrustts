use forrustts_genetic_maps::GenerateBreakpoints;
use forrustts_genetic_maps::GeneticMap;
use forrustts_genetic_maps::GeneticMapBuilder;
use rand::SeedableRng;

#[test]
fn test_build_empty() {
    let builder = GeneticMapBuilder::default();
    assert!(GeneticMap::new_from_builder(builder).is_some());
}

#[test]
fn test_build_non_empty() {
    let builder = GeneticMapBuilder::default()
        .extend_bernoullil(&[forrustts_genetic_maps::BernoulliCrossover::new(0, 1, 1.).unwrap()]);
    assert!(GeneticMap::new_from_builder(builder).is_some());
}

#[test]
fn test_generate_breakpoints() {
    let builder = GeneticMapBuilder::default()
        .extend_bernoullil(&[forrustts_genetic_maps::BernoulliCrossover::new(0, 1, 1.).unwrap()]);
    let mut map = GeneticMap::new_from_builder(builder).unwrap();
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    map.generate_breakpoints(&mut rng);
    assert_eq!(map.breakpoints().len(), 1);
}
