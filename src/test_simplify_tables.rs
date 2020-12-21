#[cfg(test)]
mod test {
    // NOTE: Currently, these tests are both testing
    // stuff from tskit_rust and forrusts, which isn't great.
    // We'll clean this up later when we get better abstractions
    // into tskit_rust.
    use crate::simplify_tables;
    use crate::tsdef::{IdType, Position, Time};
    use crate::wright_fisher::*;
    use crate::SimplificationFlags;
    use crate::SimplificationOutput;
    use crate::TableCollection;
    use std::mem::MaybeUninit;
    use tskit_rust::bindings as tskr;

    fn tables_to_treeseq(
        tables: &mut tskit_rust::TableCollection,
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

        tsk_ts
    }

    fn simulate_data(
        num_generations: Time,
        genome_length: Position,
        psurvival: f64,
        seed: usize,
        simplification_interval: Option<Time>,
        flags: SimulationFlags,
    ) -> (TableCollection, Vec<i32>) {
        // None here means "never simplify".
        neutral_wf(
            PopulationParams {
                size: 250,
                genome_length,
                littler: 5e-3,
                psurvival,
            },
            SimulationParameters {
                simplification_interval,
                seed,
                nsteps: num_generations,
                flags,
            },
        )
    }

    #[test]
    fn test_kc_distance_to_tskit() {
        // 1. Simulate in rust
        // 2. Copy to tskit and reverse time
        // 3. sort and simplify via tskit
        // 4. sort and simplify via rust
        // 5. Do we get the same stuff out?

        let num_generations = 5000;
        let genome_length = 1000000;

        let (mut tables, mut is_sample) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            42,
            None,
            SimulationFlags::empty(),
        );

        let mut tsk_tables = crate::tskit::convert_to_tskit(
            &tables,
            &is_sample,
            crate::tskit::simple_time_reverser(num_generations),
            // Do not index tables here!
            // Things are unsorted!
            false,
        );

        // Now, sort and simplify the tables we got from the sim:
        tables.sort_tables_for_simplification();
        let mut samples: Vec<IdType> = vec![];
        for (i, n) in tables.nodes().iter().enumerate() {
            if n.time == num_generations {
                samples.push(i as IdType);
            }
        }

        let mut output = SimplificationOutput::new();
        simplify_tables(
            &samples,
            SimplificationFlags::empty(),
            &mut tables,
            &mut output,
        );

        is_sample = vec![0; tables.num_nodes()];
        for i in is_sample.iter_mut().take(500) {
            *i = 1;
        }

        let mut simplified_rust_tables = crate::tskit::convert_to_tskit(
            &tables,
            &is_sample,
            crate::tskit::simple_time_reverser(num_generations),
            true,
        );

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
        }

        // Get tree sequences now
        let mut tsk_ts = tables_to_treeseq(&mut tsk_tables);
        let mut rust_ts = tables_to_treeseq(&mut simplified_rust_tables);

        unsafe {
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(tsk_ts.as_ptr()));
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
            assert!((kc - 0.).abs() < f64::EPSILON);

            tskr::tsk_treeseq_free(tsk_ts.as_mut_ptr());
            tskr::tsk_treeseq_free(rust_ts.as_mut_ptr());
        }
    }

    #[test]
    fn test_buffer_vs_sort() {
        let num_generations = 5000;
        let genome_length = 1000000;

        let flags = SimulationFlags::USE_STATE;
        let (tables_sorted, is_sample_sorted) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            14613641,
            Some(100),
            flags,
        );

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            14613641,
            Some(100),
            flags,
        );

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes_.iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes_.iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let mut tables_sorted_tskit = crate::tskit::convert_to_tskit(
            &tables_sorted,
            &is_sample_sorted,
            crate::tskit::simple_time_reverser(num_generations),
            true,
        );

        let mut tables_buffered_tskit = crate::tskit::convert_to_tskit(
            &tables_buffered,
            &is_sample_buffered,
            crate::tskit::simple_time_reverser(num_generations),
            true,
        );

        let mut sorted_ts = tables_to_treeseq(&mut tables_sorted_tskit);
        let mut buffered_ts = tables_to_treeseq(&mut tables_buffered_tskit);
        unsafe {
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(sorted_ts.as_ptr()));
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(buffered_ts.as_ptr()));
            let mut kc: f64 = -1.;
            let kcp: *mut f64 = &mut kc;
            let rv = tskr::tsk_treeseq_kc_distance(
                sorted_ts.as_mut_ptr(),
                buffered_ts.as_mut_ptr(),
                0.,
                kcp,
            );
            assert_eq!(rv, 0);
            assert!((kc - 0.).abs() < f64::EPSILON);
            tskr::tsk_treeseq_free(sorted_ts.as_mut_ptr());
            tskr::tsk_treeseq_free(buffered_ts.as_mut_ptr());
        }
    }

    // The KC distance code will barf on trees where samples
    // have different ages.  So we have to use less direct methods
    // to compare.
    #[test]
    fn test_buffer_vs_sort_overlapping_generations() {
        let num_generations = 5000;
        let genome_length = 1000000;

        let flags = SimulationFlags::USE_STATE;
        let (tables_sorted, is_sample_sorted) = simulate_data(
            num_generations,
            genome_length,
            0.5,
            14613641,
            Some(100),
            flags,
        );

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.5,
            14613641,
            Some(100),
            flags,
        );

        assert_eq!(tables_sorted.num_nodes(), tables_buffered.num_nodes());

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes_.iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes_.iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let mut tables_sorted_tskit = crate::tskit::convert_to_tskit(
            &tables_sorted,
            &is_sample_sorted,
            crate::tskit::simple_time_reverser(num_generations),
            true,
        );

        let mut tables_buffered_tskit = crate::tskit::convert_to_tskit(
            &tables_buffered,
            &is_sample_buffered,
            crate::tskit::simple_time_reverser(num_generations),
            true,
        );

        let mut sorted_ts = tables_to_treeseq(&mut tables_sorted_tskit);
        let mut buffered_ts = tables_to_treeseq(&mut tables_buffered_tskit);
        unsafe {
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(sorted_ts.as_ptr()));
            assert_eq!(500, tskr::tsk_treeseq_get_num_samples(buffered_ts.as_ptr()));
            assert_eq!(
                tskr::tsk_treeseq_get_num_trees(sorted_ts.as_ptr()),
                tskr::tsk_treeseq_get_num_trees(buffered_ts.as_ptr())
            );
            tskr::tsk_treeseq_free(sorted_ts.as_mut_ptr());
            tskr::tsk_treeseq_free(buffered_ts.as_mut_ptr());
        }
    }
}
