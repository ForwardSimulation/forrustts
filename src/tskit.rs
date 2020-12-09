//! tskit interface

use crate::tables::{TableCollection, TablesError};
use std::mem::MaybeUninit;

#[cfg(test)]
mod tests {

    use super::*;
    use tskit_rust::bindings;

    #[test]
    pub fn test_adding_rows() -> () {
        let mut tables = TableCollection::new(1000).unwrap();
        tables.add_edge(0, 1, 2, 3).unwrap();
        let mut edges: MaybeUninit<tskit_rust::bindings::tsk_edge_table_t> = MaybeUninit::uninit();
        unsafe {
            let rv = tskit_rust::bindings::tsk_edge_table_init(edges.as_mut_ptr(), 0);
            edges.assume_init();
            for e in tables.edges() {
                tskit_rust::bindings::tsk_edge_table_add_row(
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
            tskit_rust::bindings::tsk_edge_table_free(edges.as_mut_ptr());
        }
    }
}
