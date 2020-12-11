use crate::simplification_logic::{AncestryList, SegmentOverlapper};
use crate::tables::{EdgeTable, NodeTable};

/// Holds internal memory used by 
/// simplification machinery.
///
/// During simplification, several large
/// memory blocks are required. This type
/// allows those allocations to be re-used
/// in subsequent calls to
/// [simplify_tables_with_state](fn.simplify_tables_with_state.html).
/// Doing so typically improves run times at
/// the cost of higher peak memory consumption.
pub struct SimplificationBuffers {
    pub(crate) new_edges: EdgeTable,
    pub(crate) temp_edge_buffer: EdgeTable,
    pub(crate) new_nodes: NodeTable,
    pub(crate) overlapper: SegmentOverlapper,
    pub(crate) ancestry: AncestryList,
}

impl SimplificationBuffers {
    /// Create a new instance.
    pub const fn new() -> SimplificationBuffers {
        return SimplificationBuffers {
            new_edges: EdgeTable::new(),
            temp_edge_buffer: EdgeTable::new(),
            new_nodes: NodeTable::new(),
            overlapper: SegmentOverlapper::new(),
            ancestry: AncestryList::new(),
        };
    }
    
    // NOTE: should this be fully pub?
    pub(crate) fn clear(&mut self)
    {
        self.new_edges.clear();
        self.temp_edge_buffer.clear();
        self.new_nodes.clear();
    }
}
