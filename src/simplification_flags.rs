use bitflags::bitflags;

bitflags! {
    /// Boolean flags affecting simplification
    /// behavior.
    ///
    /// Currently, this is unused, and exists
    /// as a placeholder for the future.
    ///
    /// # Example
    ///
    /// ```
    /// let e = forrustts::SimplificationFlags::empty();
    /// assert!(e.contains(forrustts::SimplificationFlags::NONE));
    /// ```
    #[derive(Default)]
    pub struct SimplificationFlags: u32 {
        const NONE = 0;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_empty() {
        let e = SimplificationFlags::empty();
        assert!(e.contains(SimplificationFlags::NONE));
        assert_eq!(e.bits(), 0);
    }
}
