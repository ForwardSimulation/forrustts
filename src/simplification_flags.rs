use bitflags::bitflags;

bitflags! {
    /// Boolean flags affecting simplification
    /// behavior.
    ///
    /// # Example
    ///
    /// ```
    /// let e = forrustts::SimplificationFlags::empty();
    /// assert_eq!(e.bits(), 0);
    /// ```
    #[derive(Default)]
    pub struct SimplificationFlags: u32 {
        /// Validate that input edges are sorted
        const VALIDATE_EDGES = 1 << 0;
        /// Validate that input mutations are sorted
        const VALIDATE_MUTATIONS = 1 << 1;
        /// Validate all tables.
        const VALIDATE_ALL = Self::VALIDATE_EDGES.bits | Self::VALIDATE_MUTATIONS.bits;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_empty() {
        let e = SimplificationFlags::empty();
        assert_eq!(e.bits(), 0);
    }
}
