#![macro_use]

macro_rules! iterator_for_nodeiterator {
    ($ty: ty) => {
        impl Iterator for $ty {
            type Item = forrustts_core::newtypes::NodeId;
            fn next(&mut self) -> Option<Self::Item> {
                self.next_node();
                self.current_node()
            }
        }
    };
}
