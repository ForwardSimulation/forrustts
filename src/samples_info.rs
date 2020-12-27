use crate::IdType;

/// Information about samples used for
/// table simpilfication.
#[derive(Default)]
pub struct SamplesInfo {
    /// A list of sample IDs.
    /// Can include both "alive" and
    /// "ancient/remembered/preserved" sample
    /// nodes.
    pub samples: Vec<IdType>,
    /// When using [``EdgeBuffer``] to record transmission
    /// events, this list must contain a list of all node IDs
    /// alive the last time simplification happened. Here,
    /// "alive" means "could leave more descendants".
    /// At the *start* of a simulation, this  should be filled
    /// with a list of "founder" node IDs.
    pub edge_buffer_founder_nodes: Vec<IdType>,
}

impl SamplesInfo {
    /// Generate a new instance.
    pub fn new() -> Self {
        SamplesInfo {
            samples: vec![],
            edge_buffer_founder_nodes: vec![],
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_default() {
        let s: SamplesInfo = Default::default();
        assert!(s.samples.is_empty());
        assert!(s.edge_buffer_founder_nodes.is_empty());
    }
}
