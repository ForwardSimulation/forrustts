/// Useful information output by table
/// simplification.
pub struct SimplificationOutput {
    /// Maps input node ID to output ID.
    /// Values are set to [``NULL_ID``](crate::NULL_ID)
    /// for input nodes that "simplify out".
    pub idmap: Vec<crate::IdType>,
}

impl SimplificationOutput {
    /// Create a new instance.
    pub fn new() -> Self {
        SimplificationOutput { idmap: vec![] }
    }
}

impl Default for SimplificationOutput {
    fn default() -> Self {
        SimplificationOutput::new()
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_default() {
        let x: SimplificationOutput = Default::default();
        assert_eq!(x.idmap.is_empty(), true);
    }
}
