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

// This trait belongs in core?
// And we can give a proc macro for deriving?
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

pub trait CrossoverPosition<T>: GenomicInterval {
    fn generate_breakpoint(&self, rng: &mut T) -> Position;
}

pub trait GeneticMap<T> {
    fn generate_breakpoints(&mut self, rng: &mut T);
    fn breakpoints(&self) -> &[Position];
}

#[cfg(test)]
mod tests {
    use super::*;

    // Should be in core?
    struct Segment {
        left: Position,
        right: Position,
    }

    impl GenomicInterval for Segment {
        fn left(&self) -> Position {
            self.left
        }
        fn right(&self) -> Position {
            self.right
        }
    }

    // Should be in core?
    struct Point {
        p: Position,
    }

    impl GenomicInterval for Point {
        fn left(&self) -> Position {
            self.p
        }
        fn right(&self) -> Position {
            // See GitHub issue 240
            (i64::from(self.p) + 1).into()
        }
    }
}
