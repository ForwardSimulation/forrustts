use core::fmt::Debug;
use forrustts_core::newtypes::Position;

pub trait PoissonCrossoverRegion {
    fn mean(&self) -> f64;
}

pub trait BinomialCrossoverRegion {
    fn prob(&self) -> f64;
}

pub trait FixedNumberOfCrossoverRegion {
    fn num_breakpoints(&self) -> u32;
}

pub trait GeneticMapElementCore<T> {
    fn begin(&self) -> Position;
    fn end(&self) -> Position;
    fn generate_breakpoints(&self, rng: &mut T, breakpoints: &mut crate::PositionVec);

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

pub trait GeneticMapElement<T>: GeneticMapElementCore<T> + Send + Sync + Debug
where
    T: Send + Sync + Debug,
{
}

impl<G, T> GeneticMapElement<T> for G
where
    G: GeneticMapElementCore<T> + Send + Sync + Debug,
    T: Send + Sync + Debug,
{
}

pub trait GeneticMap<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Position];
}
