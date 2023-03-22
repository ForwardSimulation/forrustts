use forrustts_core::newtypes::Position;
use forrustts_genetic_maps::BernoulliCrossover;
use forrustts_genetic_maps::GenerateBreakpoints;
use forrustts_genetic_maps::GeneticMap;
use forrustts_genetic_maps::GeneticMapBuilder;
use forrustts_genetic_maps::IndependentAssortment;
use forrustts_genetic_maps::PoissonCrossover;
use rand::SeedableRng;

#[test]
fn test_build_empty() {
    let builder = GeneticMapBuilder::default();
    assert!(GeneticMap::new_from_builder(builder).is_some());
}

#[test]
fn test_build_non_empty() {
    let builder = GeneticMapBuilder::default()
        .extend_bernoulli(&[forrustts_genetic_maps::BernoulliCrossover::new(0, 1, 1.).unwrap()]);
    assert!(GeneticMap::new_from_builder(builder).is_some());
}

#[test]
fn test_generate_bernoulli_breakpoints() {
    let builder = GeneticMapBuilder::default()
        .extend_bernoulli(&[BernoulliCrossover::new(10, 20, 1.).unwrap()])
        .extend_bernoulli(&[BernoulliCrossover::new(0, 10, 1.).unwrap()]);
    let mut map = GeneticMap::new_from_builder(builder).unwrap();
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    map.generate_breakpoints(&mut rng);
    assert_eq!(map.breakpoints().len(), 2);
    assert!(map.breakpoints().windows(2).all(|w| { w[0] <= w[1] }));
}

#[test]
fn test_generate_poisson_breakpoints() {
    let builder = GeneticMapBuilder::default()
        .extend_poisson(&[PoissonCrossover::new(10, 20, 10.).unwrap()])
        .extend_poisson(&[PoissonCrossover::new(0, 10, 10.).unwrap()]);
    let mut map = GeneticMap::new_from_builder(builder).unwrap();
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    map.generate_breakpoints(&mut rng);
    assert!(map.breakpoints().windows(2).all(|w| { w[0] <= w[1] }));
}

#[test]
fn test_generate_poisson_breakpoints_multiple_chromosomes() {
    let builder = GeneticMapBuilder::default()
        .extend_poisson(&[PoissonCrossover::new(10, 20, 10.).unwrap()])
        .extend_poisson(&[PoissonCrossover::new(0, 10, 10.).unwrap()])
        .extend_independent_assortment(&[IndependentAssortment::new(10).unwrap()]);
    let mut map = GeneticMap::new_from_builder(builder).unwrap();
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    for _ in 0..100 {
        map.generate_breakpoints(&mut rng);
        assert!(map
            .breakpoints()
            .windows(2)
            .all(|w| { Position::from(w[0]) <= Position::from(w[1]) }));
    }
}
