use clap::{value_t, value_t_or_exit, App, Arg};
use forrustts::wright_fisher::neutral_wf;
use tskit_rust;

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
            Arg::with_name("generations")
                .short("g")
                .long("generations")
                .help("number of generations to simulate")
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
        .get_matches();

    let N = value_t_or_exit!(matches.value_of("popsize"), u32);
    let g = value_t_or_exit!(matches.value_of("generations"), i64);
    let rho = value_t_or_exit!(matches.value_of("rho"), f64);
    let simplify_input = value_t!(matches.value_of("simplification_interval"), i64).unwrap_or(-1);
    let seed = value_t_or_exit!(matches.value_of("seed"), usize);
    let outfile = value_t_or_exit!(matches.value_of("outfile"), String);

    let mut simplify: Option<i64> = None;

    if simplify_input > 0 {
        simplify = Some(simplify_input);
    }

    let r = rho / (4.0 * N as f64);

    let mut tables = neutral_wf(seed, N, g, 10000000, r, 0.0, simplify);

    let mut is_sample = vec![0 as i32; tables.num_nodes()];
    for (i, n) in tables.enumerate_nodes() {
        if n.time == g {
            is_sample[i] = 1;
        }
    }

    let mut tskit_tables = forrustts::tskit::convert_to_tskit_and_drain(
        &is_sample,
        forrustts::tskit::simple_time_reverser(g),
        if simplify.is_some() { true } else { false },
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
