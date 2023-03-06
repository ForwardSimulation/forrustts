use core::fmt::Debug;
use forrustts_core::newtypes::Position;

pub trait PoissonCrossoverRegion<T>: CrossoverPosition<T> {
    fn mean(&self) -> f64;
}

// NOTE: this is conceptually messy later, when
// we want to think about GC
pub trait BernoulliCrossoverRegion<T>: CrossoverPosition<T> {
    fn prob(&self) -> f64;
}

pub trait FixedNumberOfCrossoverRegion<T>: CrossoverPosition<T> {
    fn num_breakpoints(&self) -> u32;
}

// NOTE: Should be "single breakpoint"?
pub trait CrossoverPosition<T> {
    fn begin(&self) -> Position;
    fn end(&self) -> Position;
    fn generate_breakpoint(&self, rng: &mut T) -> Position;

    fn valid(&self, genome_length: Position) -> bool {
        self.begin() >= 0
            && self.begin() < genome_length
            && self.begin() < self.end()
            && self.end() <= genome_length
    }

    fn contains(&self, position: Position) -> bool {
        position >= self.begin() && position < self.end()
    }
}

pub trait GeneticMapElement<T>: CrossoverPosition<T> + Send + Sync + Debug
where
    T: Send + Sync + Debug,
{
}

impl<G, T> GeneticMapElement<T> for G
where
    G: CrossoverPosition<T> + Send + Sync + Debug,
    T: Send + Sync + Debug,
{
}

pub trait GeneticMap<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Position];
}
