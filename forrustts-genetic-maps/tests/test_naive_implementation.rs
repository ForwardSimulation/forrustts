use forrustts_core::newtypes::Position;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::WeightedAliasIndex;

struct PoissonInterval {
    mean: f64,
    left: Position,
    right: Position,
    uniform: rand::distributions::Uniform<i64>,
}

impl PoissonInterval {
    fn new(mean: f64, left: Position, right: Position) -> Self {
        Self {
            mean,
            left,
            right,
            uniform: rand::distributions::Uniform::new(i64::from(left), i64::from(right)),
        }
    }
}

impl<T> forrustts_genetic_maps::CrossoverPosition<T> for PoissonInterval
where
    T: Rng,
{
    fn begin(&self) -> Position {
        self.left
    }
    fn end(&self) -> Position {
        self.right
    }

    fn generate_breakpoint(&self, rng: &mut T) -> Position {
        rng.sample(self.uniform).into()
    }
}

impl<T> forrustts_genetic_maps::PoissonCrossoverRegion<T> for PoissonInterval
where
    T: Rng,
{
    fn mean(&self) -> f64 {
        self.mean
    }
}

struct PoissonGeneticMap<T>
where
    T: Rng,
{
    regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>,
    lookup: WeightedAliasIndex<f64>, // O(n) construction, O(1) lookup
    dist: rand_distr::Poisson<f64>,
    breakpoints: Vec<Position>,
}

impl<T> PoissonGeneticMap<T>
where
    T: Rng,
{
    fn new(regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>) -> Self {
        let mut weights = vec![];
        regions.iter().for_each(|i| weights.push(i.mean()));
        let total_rate = weights.iter().sum();
        let lookup = WeightedAliasIndex::new(weights).unwrap();
        let dist = rand_distr::Poisson::new(total_rate).unwrap();
        Self {
            regions,
            lookup,
            dist,
            breakpoints: vec![],
        }
    }
}

impl<T> forrustts_genetic_maps::GeneticMap<T> for PoissonGeneticMap<T>
where
    T: Rng,
{
    fn breakpoints(&self) -> &[Position] {
        &self.breakpoints
    }

    fn generate_breakpoints(&mut self, rng: &mut T) {
        self.breakpoints.clear();
        let nbreakpoints = rng.sample(self.dist) as u32;
        for _ in 0..nbreakpoints {
            let region = rng.sample(&self.lookup);
            assert!(region < self.regions.len());
            let pos = self.regions[region].generate_breakpoint(rng);
            self.breakpoints.push(pos);
        }
        self.breakpoints.sort();
    }
}

#[test]
fn test_two_regions() {
    // put the trait in scope
    use forrustts_genetic_maps::GeneticMap;
    use rand::rngs::StdRng;
    let mut rng = StdRng::seed_from_u64(42);
    let regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<StdRng>>> = vec![
        Box::new(PoissonInterval::new(10.0, 0.into(), 100.into())),
        Box::new(PoissonInterval::new(0.0, 100.into(), 200.into())),
    ];
    let mut map = PoissonGeneticMap::new(regions);

    let mut some = false;
    for _ in 0..100 {
        map.generate_breakpoints(&mut rng);
        if !map.breakpoints.is_empty() {
            some = true;
        }
        assert!(map.breakpoints.iter().all(|x| x < &100));

        // assert sorted
        assert!(map.breakpoints.windows(2).all(|w| w[0] <= w[1]));
    }

    // Here, we're just asserting that we did indeed
    // generate ANY breakpoints...
    assert!(some);
}
