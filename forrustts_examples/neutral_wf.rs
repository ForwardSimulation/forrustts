use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::wright_fisher::*;

fn main() {
    let matches = App::new("neutral_wf")
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
            Arg::with_name("rho")
                .long("rho")
                .help("Scaled rate of crossing over, 4Nr")
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
        .get_matches();

    // TODO: default params

    let popsize = value_t_or_exit!(matches.value_of("popsize"), u32);
    let g = value_t_or_exit!(matches.value_of("nsteps"), i64);
    let rho = value_t_or_exit!(matches.value_of("rho"), f64);
    let simplify_input = value_t!(matches.value_of("simplification_interval"), i64).unwrap_or(-1);
    let psurvival = value_t!(matches.value_of("psurvival"), f64).unwrap_or(0.0);
    let seed = value_t_or_exit!(matches.value_of("seed"), usize);
    let outfile = value_t_or_exit!(matches.value_of("outfile"), String);

    // TODO: parameter validation..

    let mut simplify: Option<i64> = None;

    if simplify_input > 0 {
        simplify = Some(simplify_input);
    }

    let r = rho / (4.0 * popsize as f64);

    let flags = if matches.is_present("buffer_edges") {
        SimulationFlags::BUFFER_EDGES
    } else {
        SimulationFlags::USE_STATE
    };

    let (mut tables, is_sample) = neutral_wf(
        PopulationParams {
            size: popsize,
            genome_length: 10000000,
            littler: r,
            psurvival,
        },
        SimulationParams {
            simplification_interval: simplify,
            seed,
            nsteps: g,
            flags,
        },
    )
    .unwrap();

    let mut tskit_tables = forrustts::tskit::convert_to_tskit_and_drain(
        &is_sample,
        forrustts::tskit::simple_time_reverser(g),
        simplify.is_some(),
        &mut tables,
    );

    if simplify.is_some() {
        tskit_tables.dump(&outfile, 0).unwrap();
    } else {
        tskit_tables
            .dump(&outfile, tskit_rust::TSK_NO_BUILD_INDEXES)
            .unwrap();
    }
}
