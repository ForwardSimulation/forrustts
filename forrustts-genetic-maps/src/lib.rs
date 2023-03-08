//! Building genetic maps

use forrustts_core::newtypes::Position;

pub trait PoissonCrossoverProcess<T>: CrossoverPosition<T> {
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

// This trait belongs in core
pub trait GenomicInterval {
    fn left(&self) -> Position;
    fn right(&self) -> Position;
    fn contains(&self, position: Position) -> bool {
        position >= self.left() && position < self.right()
    }
    fn valid(&self, genome_length: Position) -> bool {
        self.left() >= 0
            && self.left() < genome_length
            && self.right() < self.right()
            && self.right() <= genome_length
    }
}

pub trait CrossoverPosition<T> : GenomicInterval {
    fn generate_breakpoint(&self, rng: &mut T) -> Position;
}

pub trait GeneticMap<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Position];
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Segment {
        left: Position,
        right: Position,
    }

    struct Point {
        p: Position,
    }
}
