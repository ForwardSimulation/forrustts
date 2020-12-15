#[cfg(test)]
mod test {
    // NOTE: Currently, these tests are both testing
    // stuff from tskit_rust and forrusts, which isn't great.
    // We'll clean this up later when we get better abstractions
    // into tskit_rust.
    use crate::simplify_tables::simplify_tables;
    use crate::tables::TableCollection;
    use crate::tsdef::{SamplesVec, TsInt};
    use crate::wright_fisher::neutral_wf;
    use std::mem::MaybeUninit;
    use tskit_rust::bindings as tskr;

    fn convert_to_tskit(
        time_offset: f64,
        is_sample: &Vec<i32>,
        tables: &TableCollection,
    ) -> std::mem::MaybeUninit<tskr::tsk_table_collection_t> {
        let mut tsk_tables: MaybeUninit<tskr::tsk_table_collection_t> = MaybeUninit::uninit();
        unsafe {
            let rv = tskr::tsk_table_collection_init(tsk_tables.as_mut_ptr(), 0);
            assert_eq!(rv, 0);

            (*tsk_tables.as_mut_ptr()).sequence_length = tables.get_length() as f64;

            for (i, n) in tables.nodes().iter().enumerate() {
                let mut flag: u32 = 0;
                if is_sample[i] == 1 {
                    flag |= tskr::TSK_NODE_IS_SAMPLE;
                }
                let rv = tskr::tsk_node_table_add_row(
                    &mut (*tsk_tables.as_mut_ptr()).nodes,
                    flag,
                    -1.0 * (n.time as f64 - time_offset),
                    0,
                    tskr::TSK_NULL,
                    std::ptr::null(),
                    0,
                );
                assert_eq!(rv, i as i32);
            }

            for (i, e) in tables.edges().iter().enumerate() {
                let rv = tskr::tsk_edge_table_add_row(
                    &mut (*tsk_tables.as_mut_ptr()).edges,
                    e.left as f64,
                    e.right as f64,
                    e.parent,
                    e.child,
                    std::ptr::null(),
                    0,
                );
                assert_eq!(rv, i as i32);
            }

            // I always forget this one...
            let rv = tskr::tsk_population_table_add_row(
                &mut (*tsk_tables.as_mut_ptr()).populations,
                std::ptr::null(),
                0,
            );
            assert_eq!(rv, 0);
        }
        return tsk_tables;
    }

    fn tables_to_treeseq(
        tables: &mut MaybeUninit<tskr::tsk_table_collection_t>,
    ) -> MaybeUninit<tskr::tsk_treeseq_t> {
        let mut tsk_ts: MaybeUninit<tskr::tsk_treeseq_t> = MaybeUninit::uninit();
        unsafe {
            let rv = tskr::tsk_table_collection_build_index(tables.as_mut_ptr(), 0);
            assert_eq!(rv, 0);
            let rv = tskr::tsk_treeseq_init(
                tsk_ts.as_mut_ptr(),
                tables.as_ptr(),
                tskr::TSK_SAMPLE_LISTS,
            );
            assert_eq!(rv, 0);
        }

        return tsk_ts;
    }

    #[test]
    fn test_wright_fisher() {
        // 1. Simulate in rust
        // 2. Copy to tskit and reverse time
        // 3. sort and simplify via tskit
        // 4. sort and simplify via rust
        // 5. Do we get the same stuff out?
        let num_generations = 5000;
        let genome_length = 1000000;

        // None here means "never simplify".
        let mut tables = neutral_wf(42, 250, num_generations, genome_length, 5e-3, 0.0, None);

        let mut is_sample: Vec<i32> = vec![0; tables.num_nodes()];

        for (i, n) in tables.nodes().iter().enumerate() {
            if n.time == num_generations {
                is_sample[i] = 1;
            }
        }

        let mut tsk_tables = convert_to_tskit(num_generations as f64, &is_sample, &tables);

        // Now, sort and simplify the tables we got from the sim:
        tables.sort_tables_for_simplification();
        let mut samples = SamplesVec::new();
        for (i, n) in tables.nodes().iter().enumerate() {
            if n.time == num_generations {
                samples.push(i as TsInt);
            }
        }

        simplify_tables(&samples, &mut tables);

        is_sample = vec![0; tables.num_nodes()];
        for i in 0..500 {
            is_sample[i] = 1;
        }
        let mut simplified_rust_tables =
            convert_to_tskit(num_generations as f64, &is_sample, &tables);

        unsafe {
            let rv = tskr::tsk_table_collection_sort(tsk_tables.as_mut_ptr(), std::ptr::null(), 0);
            assert!(rv == 0);
            let rv = tskr::tsk_table_collection_simplify(
                tsk_tables.as_mut_ptr(),
                samples.as_ptr(),
                samples.len() as u32,
                0,
                std::ptr::null_mut(),
            );
            assert!(rv == 0);

            // Get tree sequences now
            let mut tsk_ts = tables_to_treeseq(&mut tsk_tables);
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(tsk_ts.as_ptr()));
            let mut rust_ts = tables_to_treeseq(&mut simplified_rust_tables);
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(rust_ts.as_ptr()));

            let ne = tskr::tsk_treeseq_get_num_edges(tsk_ts.as_ptr());
            let ne2 = tskr::tsk_treeseq_get_num_edges(rust_ts.as_ptr());
            assert_eq!(ne, ne2);
            let nn = tskr::tsk_treeseq_get_num_nodes(tsk_ts.as_ptr());
            let nn2 = tskr::tsk_treeseq_get_num_nodes(rust_ts.as_ptr());
            assert_eq!(nn, nn2);

            // We expect the rv to be 0,
            // so let's init it to something else
            let mut kc: f64 = -1.;
            let kcp: *mut f64 = &mut kc;
            let rv =
                tskr::tsk_treeseq_kc_distance(tsk_ts.as_mut_ptr(), rust_ts.as_mut_ptr(), 0., kcp);
            assert_eq!(rv, 0);
            assert_eq!(kc, 0.);

            tskr::tsk_treeseq_free(tsk_ts.as_mut_ptr());
            tskr::tsk_treeseq_free(rust_ts.as_mut_ptr());
        }
    }
}
