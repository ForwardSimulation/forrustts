#[cfg(test)]
mod test {
    use crate::simplify_tables_without_state;
    use crate::tsdef::{IdType, Position, Time};
    use crate::wright_fisher::*;
    use crate::ForrusttsError;
    use crate::SamplesInfo;
    use crate::SimplificationFlags;
    use crate::SimplificationOutput;
    use crate::TableCollection;
    use tskit::TableAccess;

    fn simulate_data(
        num_generations: Time,
        genome_length: Position,
        psurvival: f64,
        seed: usize,
        // None here means "never simplify".
        simplification_interval: Option<Time>,
        flags: SimulationFlags,
    ) -> Result<(TableCollection, Vec<i32>), ForrusttsError> {
        neutral_wf(
            PopulationParams {
                size: 250,
                genome_length,
                littler: 5e-3,
                psurvival,
            },
            SimulationParams {
                simplification_interval,
                seed,
                nsteps: num_generations,
                flags,
                simplification_flags: SimplificationFlags::VALIDATE_ALL,
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
        )
        .unwrap();

        let mut tsk_tables = crate::tskit_tools::convert_to_tskit_minimal(
            &tables,
            &is_sample,
            crate::tskit_tools::simple_time_reverser(num_generations),
            // Do not index tables here!
            // Things are unsorted!
            false,
        );

        // Now, sort and simplify the tables we got from the sim:
        tables.sort_tables(crate::TableSortingFlags::empty());
        let mut samples = SamplesInfo::new();
        for (i, n) in tables.nodes().iter().enumerate() {
            if n.time == num_generations {
                samples.samples.push(i as IdType);
            }
        }

        let mut output = SimplificationOutput::new();
        simplify_tables_without_state(
            &samples,
            SimplificationFlags::empty(),
            &mut tables,
            &mut output,
        )
        .unwrap();

        is_sample = vec![0; tables.num_nodes()];
        for i in is_sample.iter_mut().take(500) {
            *i = 1;
        }

        let simplified_rust_tables = crate::tskit_tools::convert_to_tskit_minimal(
            &tables,
            &is_sample,
            crate::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        tsk_tables
            .full_sort(tskit::TableSortOptions::default())
            .unwrap();
        tsk_tables
            .simplify(
                &samples.samples,
                tskit::SimplificationOptions::default(),
                false,
            )
            .unwrap();

        // Get tree sequences now
        let tsk_ts = tsk_tables
            .tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES)
            .unwrap();

        let rust_ts = simplified_rust_tables
            .tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES)
            .unwrap();

        assert_eq!(500, tsk_ts.num_samples());
        assert_eq!(500, rust_ts.num_samples());
        let ne = tsk_ts.edges().num_rows();
        let ne2 = rust_ts.edges().num_rows();
        assert_eq!(ne, ne2);
        let nn = tsk_ts.nodes().num_rows();
        let nn2 = rust_ts.nodes().num_rows();
        assert_eq!(nn, nn2);

        let kc = tsk_ts.kc_distance(&&rust_ts, 0.).unwrap();
        assert!((kc - 0.).abs() < f64::EPSILON);
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
        )
        .unwrap();

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.0,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes_.iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes_.iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let tables_sorted_tskit = crate::tskit_tools::convert_to_tskit_minimal(
            &tables_sorted,
            &is_sample_sorted,
            crate::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let tables_buffered_tskit = crate::tskit_tools::convert_to_tskit_minimal(
            &tables_buffered,
            &is_sample_buffered,
            crate::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let sorted_ts = tables_sorted_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        let buffered_ts = tables_buffered_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        assert_eq!(500, sorted_ts.num_samples());
        assert_eq!(500, buffered_ts.num_samples());
        let kc = sorted_ts.kc_distance(&&buffered_ts, 0.).unwrap();
        assert!((kc - 0.).abs() < f64::EPSILON);
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
        )
        .unwrap();

        let flags = SimulationFlags::BUFFER_EDGES;
        let (tables_buffered, is_sample_buffered) = simulate_data(
            num_generations,
            genome_length,
            0.5,
            14613641,
            Some(100),
            flags,
        )
        .unwrap();

        assert_eq!(tables_sorted.num_nodes(), tables_buffered.num_nodes());

        // The sums of node times should be the same if we've got
        // both methods working
        let sum_times_sorted: Time = tables_sorted.nodes_.iter().map(|x| x.time).sum();
        let sum_times_buffered: Time = tables_buffered.nodes_.iter().map(|x| x.time).sum();
        assert_eq!(sum_times_sorted, sum_times_buffered);

        let tables_sorted_tskit = crate::tskit_tools::convert_to_tskit_minimal(
            &tables_sorted,
            &is_sample_sorted,
            crate::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let tables_buffered_tskit = crate::tskit_tools::convert_to_tskit_minimal(
            &tables_buffered,
            &is_sample_buffered,
            crate::tskit_tools::simple_time_reverser(num_generations),
            true,
        );

        let sorted_ts = tables_sorted_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        let buffered_ts = tables_buffered_tskit
            .tree_sequence(tskit::TreeSequenceFlags::default())
            .unwrap();
        assert_eq!(500, sorted_ts.num_samples());
        assert_eq!(500, buffered_ts.num_samples());
        assert_eq!(sorted_ts.num_trees(), buffered_ts.num_trees());
    }
}
