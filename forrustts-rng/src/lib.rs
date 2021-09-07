//! Random number generation.
//!
//! This crate is a set of thin abstractions
//! around the [`rgsl`](https://docs.rs/GSL/) crate.

/// A random number generator.
///
/// This is a newtype wrapper
/// around the GNU Scientific Library
/// mersenne twister (mt19937).
#[repr(transparent)]
pub struct Rng(rgsl::Rng);

impl Rng {
    /// Create a new [`Rng`] with a seed.
    pub fn new(seed: usize) -> Self {
        let mut rng = rgsl::rng::Rng::new(rgsl::rng::algorithms::mt19937()).unwrap();
        rng.set(seed);

        Self(rng)
    }
}

/// Provide access to the underlying rng type
/// wrapped by [`Rng`].
pub trait UnderlyingRngAccess {
    type UnderlyingRng;
    /// Get a reference to the underlying rng
    fn as_underlying_ref(&self) -> &Self::UnderlyingRng;
    /// Get a mutable reference to the underlying rng
    fn as_underlying_mut_ref(&mut self) -> &mut Self::UnderlyingRng;
}

impl UnderlyingRngAccess for Rng {
    type UnderlyingRng = rgsl::rng::Rng;
    fn as_underlying_ref(&self) -> &Self::UnderlyingRng {
        &self.0
    }
    fn as_underlying_mut_ref(&mut self) -> &mut Self::UnderlyingRng {
        &mut self.0
    }
}

/// Return a Poisson deviate from a distribution
/// with a given `mean`.
///
/// # Example
///
/// ```
/// let mut rng = forrustts_rng::Rng::new(42);
/// let _ = forrustts_rng::poisson(&mut rng, 5e-3);
/// ```
#[inline]
pub fn poisson(rng: &mut Rng, mean: f64) -> u32 {
    rng.0.poisson(mean)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mut_ref() {
        use UnderlyingRngAccess;
        let mut rng = Rng::new(101);
        let _ = rng.as_underlying_mut_ref().poisson(10.);
    }
}
