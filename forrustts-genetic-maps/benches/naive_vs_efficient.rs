#[macro_use]
extern crate bencher;

use bencher::Bencher;

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

struct PoissonMapBuilder<T>
where
    T: Rng,
{
    regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>,
}

impl<T> Default for PoissonMapBuilder<T>
where
    T: Rng,
{
    fn default() -> Self {
        Self { regions: vec![] }
    }
}

impl<T> PoissonMapBuilder<T>
where
    T: Rng,
{
    fn add_region<P>(&mut self, region: P)
    where
        P: forrustts_genetic_maps::PoissonCrossoverRegion<T> + 'static,
    {
        self.regions.push(Box::new(region));
    }

    fn build(self) -> PoissonGeneticMap<T> {
        PoissonGeneticMap::new(self.regions)
    }
}

struct NaivePoissonGeneticMap<T>
where
    T: Rng,
{
    regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>,
    poissons: Vec<rand_distr::Poisson<f64>>,
    breakpoints: Vec<Position>,
}

impl<T> NaivePoissonGeneticMap<T>
where
    T: Rng,
{
    fn new(regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>) -> Self {
        let mut poissons = vec![];
        for r in regions.iter() {
            if r.mean() > 0.0 {
                let p = rand_distr::Poisson::new(r.mean()).unwrap();
                poissons.push(p);
            }
        }
        Self {
            regions,
            poissons,
            breakpoints: vec![],
        }
    }
}

impl<T> forrustts_genetic_maps::GeneticMap<T> for NaivePoissonGeneticMap<T>
where
    T: Rng,
{
    fn breakpoints(&self) -> &[Position] {
        &self.breakpoints
    }

    fn generate_breakpoints(&mut self, rng: &mut T) {
        self.breakpoints.clear();
        for (i, p) in self.poissons.iter().enumerate() {
            let n = rng.sample(p) as u32;
            for _ in 0..n {
                self.breakpoints
                    .push(self.regions[i].generate_breakpoint(rng));
            }
        }
        self.breakpoints.sort();
    }
}

struct NaivePoissonMapBuilder<T>
where
    T: Rng,
{
    regions: Vec<Box<dyn forrustts_genetic_maps::PoissonCrossoverRegion<T>>>,
}

impl<T> Default for NaivePoissonMapBuilder<T>
where
    T: Rng,
{
    fn default() -> Self {
        Self { regions: vec![] }
    }
}

impl<T> NaivePoissonMapBuilder<T>
where
    T: Rng,
{
    fn add_region<P>(&mut self, region: P)
    where
        P: forrustts_genetic_maps::PoissonCrossoverRegion<T> + 'static,
    {
        self.regions.push(Box::new(region));
    }

    fn build(self) -> NaivePoissonGeneticMap<T> {
        NaivePoissonGeneticMap::new(self.regions)
    }
}

fn efficient(total_rate: f64, num_regions: u32) {
    // put the trait in scope
    use forrustts_genetic_maps::GeneticMap;
    use rand::rngs::StdRng;
    let mut rng = StdRng::seed_from_u64(42);

    let mut builder = PoissonMapBuilder::<StdRng>::default();
    let rate_per_region = total_rate / (num_regions as f64);
    for _ in 0..num_regions {
        builder.add_region(PoissonInterval::new(rate_per_region, 0.into(), 100.into()));
    }
    let mut map = builder.build();

    for _ in 0..100 {
        map.generate_breakpoints(&mut rng);
    }
}

fn naive(total_rate: f64, num_regions: u32) {
    // put the trait in scope
    use forrustts_genetic_maps::GeneticMap;
    use rand::rngs::StdRng;
    let mut rng = StdRng::seed_from_u64(42);

    let mut builder = NaivePoissonMapBuilder::<StdRng>::default();
    let rate_per_region = total_rate / (num_regions as f64);
    for _ in 0..num_regions {
        builder.add_region(PoissonInterval::new(rate_per_region, 0.into(), 100.into()));
    }
    let mut map = builder.build();

    for _ in 0..100 {
        map.generate_breakpoints(&mut rng);
        if !map.breakpoints.is_empty() {}
    }
}

fn run_efficient_100(bench: &mut Bencher) {
    bench.iter(|| efficient(10., 100));
}

fn run_naive_100(bench: &mut Bencher) {
    bench.iter(|| naive(10., 100));
}

fn run_efficient_1000(bench: &mut Bencher) {
    bench.iter(|| efficient(10., 1000));
}

fn run_naive_1000(bench: &mut Bencher) {
    bench.iter(|| naive(10., 1000));
}

benchmark_group!(
    benches,
    run_efficient_100,
    run_naive_100,
    run_naive_1000,
    run_efficient_1000
);
benchmark_main!(benches);
