use forrustts_tskit::simple_time_reverser;

#[derive(serde::Serialize, serde::Deserialize, Eq, PartialEq, Copy, Clone, Debug)]
struct Metadata {
    a: i32,
    b: i64,
}

impl tskit::metadata::MetadataRoundtrip for Metadata {
    fn encode(&self) -> Result<Vec<u8>, tskit::metadata::MetadataError> {
        match serde_json::to_string(self) {
            Ok(x) => Ok(x.as_bytes().to_vec()),
            Err(e) => Err(tskit::metadata::MetadataError::RoundtripError { value: Box::new(e) }),
        }
    }

    fn decode(md: &[u8]) -> Result<Self, tskit::metadata::MetadataError> {
        let value: Result<Self, serde_json::Error> = serde_json::from_slice(md);
        match value {
            Ok(v) => Ok(v),
            Err(e) => Err(tskit::metadata::MetadataError::RoundtripError { value: Box::new(e) }),
        }
    }
}

impl tskit::metadata::MutationMetadata for Metadata {}
impl tskit::metadata::SiteMetadata for Metadata {}
impl tskit::metadata::PopulationMetadata for Metadata {}
impl tskit::metadata::NodeMetadata for Metadata {}

struct TablesWithMetadata {
    tables: forrustts_tables_trees::TableCollection,
    metadata: Vec<Metadata>,
    metadata_hash: std::collections::HashMap<usize, Metadata>,
    metadata_fnvhash: fnv::FnvHashMap<usize, Metadata>,
}

fn make_reverser() -> impl Fn(forrustts_tables_trees::Time) -> f64 {
    simple_time_reverser(1)
}

impl Default for TablesWithMetadata {
    fn default() -> Self {
        let mut tables = forrustts_tables_trees::TableCollection::new(10000).unwrap();

        let node_id = tables.add_node(1, 0).unwrap();
        let site_id = tables.add_site(500, None).unwrap();
        let mutation_id = tables
            .add_mutation(node_id, Some(0), site_id, 0, None, true)
            .unwrap();

        let metadata = vec![Metadata { a: -3, b: 99 }];
        let mut metadata_hash = std::collections::HashMap::new();
        metadata_hash.insert(usize::from(mutation_id), metadata[usize::from(mutation_id)]);
        let mut metadata_fnvhash = fnv::FnvHashMap::default();
        metadata_fnvhash.insert(usize::from(mutation_id), metadata[usize::from(mutation_id)]);

        assert_eq!(
            metadata[usize::from(mutation_id)],
            metadata_hash[&usize::from(mutation_id)]
        );

        Self {
            tables,
            metadata,
            metadata_hash,
            metadata_fnvhash,
        }
    }
}

macro_rules! build_metadata_roundtrip_test {
    ($testfn: ident, $input_table: ident, $output_table: ident, $export_fn: ident) => {
        #[test]
        fn $testfn() {
            use forrustts_tables_trees::TableTypeIntoRaw;
            use tskit::TableAccess;
            let data = TablesWithMetadata::default();

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(data, $input_table, tsk_tables, metadata, $export_fn);
                validate_metadata!(tsk_tables, $output_table, data, metadata, 0, 0);
            }

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(
                    data,
                    $input_table,
                    tsk_tables,
                    metadata_hash,
                    $export_fn
                );
                validate_metadata!(tsk_tables, $output_table, data, metadata_hash, 0, &0);
            }

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(
                    data,
                    $input_table,
                    tsk_tables,
                    metadata_fnvhash,
                    $export_fn
                );
                validate_metadata!(tsk_tables, $output_table, data, metadata_fnvhash, 0, &0);
            }
        }
    };

    ($testfn: ident, $input_table: ident, $output_table: ident, $export_fn: ident, $callback: ident) => {
        #[test]
        fn $testfn() {
            use forrustts_tables_trees::TableTypeIntoRaw;
            use tskit::TableAccess;
            let data = TablesWithMetadata::default();

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(
                    data,
                    $input_table,
                    tsk_tables,
                    metadata,
                    $export_fn,
                    $callback
                );
                validate_metadata!(tsk_tables, $output_table, data, metadata, 0, 0);
            }

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(
                    data,
                    $input_table,
                    tsk_tables,
                    metadata_hash,
                    $export_fn,
                    $callback
                );
                validate_metadata!(tsk_tables, $output_table, data, metadata_hash, 0, &0);
            }

            {
                let mut tsk_tables =
                    tskit::TableCollection::new(data.tables.genome_length().into_raw() as f64)
                        .unwrap();
                export_table_with_metadata!(
                    data,
                    $input_table,
                    tsk_tables,
                    metadata_fnvhash,
                    $export_fn,
                    $callback
                );
                validate_metadata!(tsk_tables, $output_table, data, metadata_fnvhash, 0, &0);
            }
        }
    };
}

macro_rules! export_table_with_metadata {
    ($object: ident, $table: ident, $tsk_tables: ident, $metadata: ident, $fn: ident) => {
        forrustts_tskit::$fn(
            $object.tables.$table(),
            &$object.$metadata,
            &mut $tsk_tables,
        )
        .unwrap();
    };

    ($object: ident, $table: ident, $tsk_tables: ident, $metadata: ident, $fn: ident, $callback: ident) => {
        forrustts_tskit::$fn(
            $object.tables.$table(),
            &$callback(),
            &$object.$metadata,
            &mut $tsk_tables,
        )
        .unwrap();
    };
}

macro_rules! validate_metadata {
    ($tsk_tables: ident, $table: ident, $object: ident, $container: ident, $index: expr, $index_two: expr) => {
        match $tsk_tables
            .$table()
            .metadata::<Metadata>($index.into())
            .unwrap()
        {
            Some(x) => assert_eq!(x, $object.$container[$index_two]),
            None => panic!("expeced Some(metadata)"),
        }
    };
}

build_metadata_roundtrip_test!(
    test_mutation_metadata_roundtrips,
    mutations,
    mutations,
    export_mutations_with_metadata,
    make_reverser
);

build_metadata_roundtrip_test!(
    test_site_metadata_roundtrips,
    sites,
    sites,
    export_sites_with_metadata
);

build_metadata_roundtrip_test!(
    test_population_metadata_roundtrips,
    nodes,
    populations,
    build_population_table_with_metadata
);

build_metadata_roundtrip_test!(
    test_node_metadata_roundtrips,
    nodes,
    nodes,
    export_nodes_with_metadata,
    make_reverser
);