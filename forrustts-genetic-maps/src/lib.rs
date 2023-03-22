//! Core types for genetic maps

use forrustts_core::newtypes::Position;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum Breakpoint {
    Crossover(Position),
    IndependentAssortment(Position),
}

pub trait GenerateBreakpoints<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Breakpoint];
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PoissonCrossover {
    left: Position,
    right: Position,
    mean: f64,
}

impl PoissonCrossover {
    pub fn new<L, R>(left: L, right: R, mean: f64) -> Option<Self>
    where
        L: Into<Position>,
        R: Into<Position>,
    {
        let left = left.into();
        let right = right.into();
        if left < 0 || right < 0 || right <= left || !mean.is_finite() || mean < 0.0 {
            None
        } else {
            Some(Self { left, right, mean })
        }
    }

    pub fn left(&self) -> Position {
        self.left
    }
    pub fn right(&self) -> Position {
        self.right
    }
    pub fn mean(&self) -> f64 {
        self.mean
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BernoulliCrossover {
    left: Position,
    right: Position,
    probability: f64,
}

impl BernoulliCrossover {
    pub fn new<L, R>(left: L, right: R, probability: f64) -> Option<Self>
    where
        L: Into<Position>,
        R: Into<Position>,
    {
        let left = left.into();
        let right = right.into();
        if left < 0
            || right < 0
            || right <= left
            || !probability.is_finite()
            || !(0.0..=1.0).contains(&probability)
        {
            None
        } else {
            Some(Self {
                left,
                right,
                probability,
            })
        }
    }
    pub fn left(&self) -> Position {
        self.left
    }
    pub fn right(&self) -> Position {
        self.right
    }
    pub fn probability(&self) -> f64 {
        self.probability
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IndependentAssortment {
    at: Position,
}

impl IndependentAssortment {
    pub fn new<P: Into<Position>>(at: P) -> Option<Self> {
        let at = at.into();
        if at >= 0 {
            Some(Self { at })
        } else {
            None
        }
    }

    pub fn at(&self) -> Breakpoint {
        Breakpoint::IndependentAssortment(self.at)
    }
}

#[derive(Default, Debug, Clone)]
pub struct GeneticMapBuilder {
    pub poisson: Vec<PoissonCrossover>,
    pub binomial: Vec<BernoulliCrossover>,
    pub independent_assortment: Vec<IndependentAssortment>,
}

impl GeneticMapBuilder {
    pub fn extend_poisson(mut self, intervals: &[PoissonCrossover]) -> Self {
        self.poisson.extend_from_slice(intervals);
        self
    }
    pub fn extend_binomial(mut self, intervals: &[BernoulliCrossover]) -> Self {
        self.binomial.extend_from_slice(intervals);
        self
    }
    pub fn extend_independent_assortment(mut self, at: &[IndependentAssortment]) -> Self {
        self.independent_assortment.extend_from_slice(at);
        self
    }

    pub fn poisson(&self) -> &[PoissonCrossover] {
        &self.poisson
    }

    pub fn binomial(&self) -> &[BernoulliCrossover] {
        &self.binomial
    }

    pub fn independent_assortment(&self) -> &[IndependentAssortment] {
        &self.independent_assortment
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyGeneticMap {}

    impl GenerateBreakpoints<()> for MyGeneticMap {
        fn generate_breakpoints(&mut self, _rng: &mut ()) {}
        fn breakpoints(&self) -> &[Breakpoint] {
            &[]
        }
    }
}

#[test]
fn test_poisson_crossover() {
    assert!(PoissonCrossover::new(0, 1, 1e-3).is_some());
    assert!(PoissonCrossover::new(0, 1, -1e-3).is_none());
    assert!(PoissonCrossover::new(0, 1, f64::NEG_INFINITY).is_none());
    assert!(PoissonCrossover::new(0, 1, f64::NAN).is_none());
    assert!(PoissonCrossover::new(1, 0, 1e-3).is_none());
}

#[test]
fn test_binomial_crossover() {
    assert!(BernoulliCrossover::new(0, 1, 1e-3).is_some());
    assert!(BernoulliCrossover::new(0, 1, 1.0).is_some());
    assert!(BernoulliCrossover::new(0, 1, 1.0 + f64::EPSILON).is_none());
    assert!(BernoulliCrossover::new(0, 1, -1e-3).is_none());
    assert!(BernoulliCrossover::new(0, 1, f64::NEG_INFINITY).is_none());
    assert!(BernoulliCrossover::new(0, 1, f64::NAN).is_none());
    assert!(BernoulliCrossover::new(1, 0, 1e-3).is_none());
}

#[test]
fn test_independent_assortment() {
    assert!(IndependentAssortment::new(1).is_some());
    assert!(IndependentAssortment::new(0).is_some());
    assert!(IndependentAssortment::new(-1).is_none());
}

#[test]
fn test_builder() {
    let builder = GeneticMapBuilder::default()
        .extend_poisson(&[PoissonCrossover::new(0, 1, 1e-3).unwrap()])
        .extend_binomial(&[BernoulliCrossover::new(0, 1, 0.25).unwrap()])
        .extend_poisson(&[PoissonCrossover::new(1, 2, 2e-3).unwrap()])
        .extend_independent_assortment(&[IndependentAssortment::new(1).unwrap()]);
    assert_eq!(builder.poisson().len(), 2);
    assert_eq!(builder.binomial().len(), 1);
    assert_eq!(builder.independent_assortment().len(), 1);
}
