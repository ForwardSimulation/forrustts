//! Building genetic maps

use forrustts_core::newtypes::Position;
use rand::Rng;

use core::fmt::Debug;

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

/// A newtype wrapper around `Vec<Position>`
/// that only provides `append`.
#[derive(Debug, Default)]
#[repr(transparent)]
pub struct PositionVec(Vec<Position>);

impl PositionVec {
    pub fn push(&mut self, position: Position) {
        self.0.push(position);
    }
}

#[derive(Debug)]
pub struct BoxedGeneticMap<T>
where
    T: Send + Sync,
{
    map: Vec<Box<dyn GeneticMapElement<T>>>,
    breakpoints: PositionVec,
}

impl<T> BoxedGeneticMap<T>
where
    T: Send + Sync,
{
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_from_vec(map: Vec<Box<dyn GeneticMapElement<T>>>) -> Self {
        Self {
            map,
            breakpoints: PositionVec::default(),
        }
    }

    pub fn add_element(&mut self, element: Box<dyn GeneticMapElement<T>>) {
        self.map.push(element);
    }
}

impl<T> Default for BoxedGeneticMap<T>
where
    T: Send + Sync,
{
    fn default() -> Self {
        Self {
            map: vec![],
            breakpoints: PositionVec::default(),
        }
    }
}

impl<T> GeneticMap<T> for BoxedGeneticMap<T>
where
    T: Send + Sync,
{
    fn generate_breakpoints(&mut self, rng: &mut T) {
        self.breakpoints.0.clear();
        for i in &self.map {
            i.generate_breakpoints(rng, &mut self.breakpoints);
        }
        self.breakpoints.0.sort();
        // TODO: eliminate need for sentinel
        if !self.breakpoints.0.is_empty() {
            self.breakpoints.push(Position::MAX);
        }
    }

    fn breakpoints(&self) -> &[Position] {
        &self.breakpoints.0
    }
}

#[derive(Copy, Clone, Debug)]
pub struct PoissonInterval {
    beg: Position,
    end: Position,
    mean: f64,
    uniform: rand::distributions::Uniform<i64>,
}

impl PoissonInterval {
    pub fn new<P: Into<Position> + Copy>(beg: P, end: P, mean: f64) -> Option<Self> {
        if mean < 0.0 || !mean.is_finite() {
            None
        } else {
            Some(Self {
                beg: beg.into(),
                end: end.into(),
                mean,
                uniform: rand::distributions::Uniform::new(
                    i64::from(beg.into()),
                    i64::from(end.into()),
                ),
            })
        }
    }
}

impl std::fmt::Display for PoissonInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PoissonInterval(beg = {}, end = {}, mean = {})",
            i64::from(self.beg),
            i64::from(self.end),
            self.mean
        )
    }
}

impl<T> CrossoverPosition<T> for PoissonInterval
where
    T: Rng,
{
    fn begin(&self) -> Position {
        self.beg
    }

    fn end(&self) -> Position {
        self.end
    }

    fn generate_breakpoint(&self, rng: &mut T) -> Position {
        rng.sample(self.uniform).into()
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BinomialPoint {
    position: Position,
    probability: f64,
    dist: rand::distributions::Bernoulli,
}

impl BinomialPoint {
    pub fn new<P: Into<Position>>(position: P, probability: f64) -> Option<Self> {
        Some(Self {
            position: position.into(),
            probability,
            // FIXME
            dist: rand::distributions::Bernoulli::new(probability).ok()?,
        })
    }
}

impl std::fmt::Display for BinomialPoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "BinomialPoint(position = {}, probability = {})",
            i64::from(self.position),
            self.probability,
        )
    }
}

impl<T> CrossoverPosition<T> for BinomialPoint
where
    T: Rng,
{
    fn begin(&self) -> Position {
        self.position
    }

    // FIXME: probably wrong
    fn end(&self) -> Position {
        self.position
    }

    fn generate_breakpoint(&self, _rng: &mut T) -> Position {
        self.position
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BinomialInterval {
    beg: Position,
    end: Position,
    probability: f64,
    dist: rand::distributions::Bernoulli,
    uniform: rand::distributions::Uniform<i64>,
}

impl BinomialInterval {
    pub fn new<P: Into<Position> + Copy>(beg: P, end: P, probability: f64) -> Option<Self> {
        Some(Self {
            beg: beg.into(),
            end: end.into(),
            probability,
            dist: rand::distributions::Bernoulli::new(probability).ok()?,
            uniform: rand::distributions::Uniform::new(
                i64::from(beg.into()),
                i64::from(end.into()),
            ),
        })
    }
}

impl std::fmt::Display for BinomialInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "BinomialPoint(beg = {}, end = {}, probability = {})",
            i64::from(self.beg),
            i64::from(self.end),
            self.probability,
        )
    }
}

impl<T> CrossoverPosition<T> for BinomialInterval
where
    T: Rng,
{
    fn begin(&self) -> Position {
        self.beg
    }

    fn end(&self) -> Position {
        self.end
    }

    fn generate_breakpoint(&self, rng: &mut T) -> Position {
        rng.sample(self.uniform).into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn is_normal<T: Send + Sync>(_: T) {}

    #[test]
    fn test_types_are_normal() {
        let rng = rand::rngs::StdRng::seed_from_u64(0);
        is_normal(rng);
        is_normal(BinomialPoint::new(1, 0.5).unwrap());
        is_normal(BinomialInterval::new(0, 100, 0.3).unwrap());
        let p = PoissonInterval::new(0, 100, 1e-3).unwrap();
        is_normal(p);

        let v: Vec<Box<dyn GeneticMapElement<rand::rngs::StdRng>>> = vec![Box::new(p)];

        let m = BoxedGeneticMap::new_from_vec(v);
        is_normal(m);
    }
}
