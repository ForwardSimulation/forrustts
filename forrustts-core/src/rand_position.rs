use rand::distributions::uniform::{SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::prelude::Rng;

use crate::Position;

#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
pub struct UniformPos(UniformInt<i64>);

impl UniformSampler for UniformPos {
    type X = Position;
    fn new<B1, B2>(low: B1, high: B2) -> Self
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformPos(UniformInt::<i64>::new(
            i64::from(*low.borrow()),
            i64::from(*high.borrow()),
        ))
    }
    fn new_inclusive<B1, B2>(low: B1, high: B2) -> Self
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformPos(UniformInt::<i64>::new_inclusive(
            i64::from(*low.borrow()),
            i64::from(*high.borrow()),
        ))
    }
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        Position::new_valid(self.0.sample(rng))
    }
}

impl SampleUniform for Position {
    type Sampler = UniformPos;
}
