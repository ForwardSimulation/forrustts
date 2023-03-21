//! Core types for genetic maps

use forrustts_core::newtypes::Position;

#[non_exhaustive]
pub enum Breakpoint {
    Crossover(Position),
    IndependentAssortment(Position),
}

pub trait GeneticMap<T> {
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
    pub fn new(left: Position, right: Position, mean: f64) -> Option<Self> {
        if left < 0 || right < 0 || right <= left || !mean.is_finite() || mean < 0.0 {
            None
        } else {
            Some(Self { left, right, mean })
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BinomialCrossover {
    left: Position,
    right: Position,
    probability: f64,
}

impl BinomialCrossover {
    pub fn new(left: Position, right: Position, probability: f64) -> Option<Self> {
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
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct IndependentAssortment {
    at: Position,
}

impl IndependentAssortment {
    pub fn new(at: Position) -> Option<Self> {
        if at >= 0 {
            Some(Self { at })
        } else {
            None
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct GeneticMapBuilder {
    poisson: Vec<PoissonCrossover>,
    binomial: Vec<BinomialCrossover>,
    independent_assortment: Vec<IndependentAssortment>,
}

impl GeneticMapBuilder {
    pub fn extend_poisson(mut self, intervals: &[PoissonCrossover]) -> Self {
        self.poisson.extend_from_slice(intervals);
        self
    }
    pub fn extend_binomial(mut self, intervals: &[BinomialCrossover]) -> Self {
        self.binomial.extend_from_slice(intervals);
        self
    }
    pub fn extend_independent_assortment(mut self, at: &[IndependentAssortment]) -> Self {
        self.independent_assortment.extend_from_slice(at);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MyGeneticMap {}

    impl GeneticMap<()> for MyGeneticMap {
        fn generate_breakpoints(&mut self, _rng: &mut ()) {}
        fn breakpoints(&self) -> &[Breakpoint] {
            &[]
        }
    }
}
