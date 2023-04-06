#[cfg(feature = "rand")]
mod test_rand_traits {
    use forrustts_core::Position;
    use proptest::prelude::*;
    use rand::SeedableRng;

    proptest! {
        #[test]
        fn test_uniform_position(a in 0..i64::MAX, b in 0..i64::MAX,
                                 seed in 0..u64::MAX) {
            if a != b { // else rand will panic
                let lo = Position::new_valid(std::cmp::min(a, b));
                let hi = Position::new_valid(std::cmp::max(a, b));
                let upos = rand::distributions::Uniform::<Position>::new(lo, hi);
                let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
                for _ in 0..100 {
                    let _ = rng.sample(upos);
                }
            }
        }
    }

    proptest! {
        #[test]
        fn test_uniform_position_inclusive(a in 0..i64::MAX, b in 0..i64::MAX,
                                           seed in 0..u64::MAX) {
            if a != b { // else rand will panic
                let lo = Position::new_valid(std::cmp::min(a, b));
                let hi = Position::new_valid(std::cmp::max(a, b));
                let upos = rand::distributions::Uniform::<Position>::new_inclusive(lo, hi);
                let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
                for _ in 0..100 {
                    let _ = rng.sample(upos);
                }
            }
        }
    }
}
