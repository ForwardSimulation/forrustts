//! tskit interface

#[cfg(test)]
mod tests {

    use crate::tables::TableCollection;
    use std::mem::MaybeUninit;
    use tskit_rust::bindings as tskr;

    #[test]
    pub fn test_adding_rows() -> () {
        let mut tables = TableCollection::new(1000).unwrap();
        tables.add_edge(0, 1, 2, 3).unwrap();
        let mut edges: MaybeUninit<tskr::tsk_edge_table_t> = MaybeUninit::uninit();
        unsafe {
            let rv = tskr::tsk_edge_table_init(edges.as_mut_ptr(), 0);
            assert_eq!(rv, 0);
            edges.assume_init();
            for e in tables.edges() {
                tskr::tsk_edge_table_add_row(
                    edges.as_mut_ptr(),
                    e.left as f64,
                    e.right as f64,
                    e.parent,
                    e.child,
                    std::ptr::null(),
                    0,
                );
            }
            assert_eq!((*edges.as_ptr()).num_rows, 1);
            tskr::tsk_edge_table_free(edges.as_mut_ptr());
        }
    }
}
