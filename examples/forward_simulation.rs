#[path = "../tests/neutral_wright_fisher.rs"]
mod neutral_wright_fisher;

use clap::{value_t, value_t_or_exit, App, Arg};
use neutral_wright_fisher::*;

fn main() {
    let matches = App::new("forward_simulation")
        .arg(
            Arg::with_name("popsize")
                .short("N")
                .long("popsize")
                .help("Diploid population size")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("nsteps")
                .short("n")
                .long("nsteps")
                .help("number of birth steps to simulate")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("xovers")
                .short("x")
                .long("xovers")
                .help("Mean number of crossovers per meiosis.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mutrate")
                .short("m")
                .long("mutrate")
                .help("Mean number of mutations per new gamete.")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("simplification_interval")
                .short("s")
                .long("simplify")
                .help("number of generations between simplifications")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("outfile")
                .short("o")
                .long("outfile")
                .help("Name of output file. The format is a tskit \"trees\" file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("seed")
                .short("S")
                .long("seed")
                .help("Random number seed")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("buffer_edges")
                .short("B")
                .long("buffer_edges")
                .help("Use edge buffering instead of sorting")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("psurvival")
                .short("P")
                .long("psurvival")
                .help("Survival probability")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("validate_tables")
                .short("v")
                .long("validate_tables")
                .help("Validate all tables prior to simplification")
                .takes_value(false),
        )
        .get_matches();

    // TODO: default params

    let popsize = value_t_or_exit!(matches.value_of("popsize"), u32);
    let nsteps = value_t_or_exit!(matches.value_of("nsteps"), i64);
    let xovers = value_t_or_exit!(matches.value_of("xovers"), f64);
    let mutrate = value_t_or_exit!(matches.value_of("mutrate"), f64);
    let simplify_input = value_t!(matches.value_of("simplification_interval"), i64).unwrap_or(-1);
    let psurvival = value_t!(matches.value_of("psurvival"), f64).unwrap_or(0.0);
    let seed = value_t_or_exit!(matches.value_of("seed"), u64);
    let outfile = value_t_or_exit!(matches.value_of("outfile"), String);
    let validate_tables = matches.is_present("validate_tables");

    // TODO: parameter validation..

    let mut simplify: Option<i64> = None;

    if simplify_input > 0 {
        simplify = Some(simplify_input);
    }

    let mut simplification_flags = forrustts::SimplificationFlags::empty();

    if validate_tables {
        simplification_flags |= forrustts::SimplificationFlags::VALIDATE_ALL;
    }

    let simparams = SimulationParams {
        popsize,
        mutrate,
        psurvival,
        xovers,
        genome_length: 10000000.into(),
        buffer_edges: true,
        simplification_interval: Some(100),
        seed,
        nsteps,
        flags: SimulationFlags::USE_STATE | SimulationFlags::BUFFER_EDGES,
        simplification_flags: forrustts::SimplificationFlags::empty(),
    };

    let (mut tables, is_sample) = neutral_wf(simparams).unwrap();

    let mut tskit_tables = forrustts::tskit_tools::convert_to_tskit_and_drain_minimal(
        &is_sample,
        forrustts::tskit_tools::simple_time_reverser(nsteps.into()),
        simplify.is_some(),
        &mut tables,
    );

    add_tskit_mutation_site_tables(&tables, nsteps.into(), &mut tskit_tables);

    let ts = match tskit_tables.tree_sequence(tskit::TreeSequenceFlags::BUILD_INDEXES) {
        Ok(x) => x,
        Err(e) => panic!("{}", e),
    };

    ts.dump(&outfile, tskit::TableOutputOptions::default())
        .unwrap();
}
