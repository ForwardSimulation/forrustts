#!/bin/bash

cargo build --all-targets
cargo build --all-targets --release

./target/release/examples/forward_simulation --popsize 5000 --nsteps 5000 -x 1e-3 -S 10101 -s 100 -o regular.trees --mutrate 0. > /dev/null
./target/release/examples/forward_simulation --popsize 5000 --nsteps 5000 -x 1e-3 -S 10101 -s 100 -o consume.trees --mutrate 0. --threads > /dev/null
# ./target/debug/examples/forward_simulation --popsize 2 --nsteps 2 -x 0. -S 101 -s 1 -o regular.trees --mutrate 0.
# ./target/debug/examples/forward_simulation --popsize 2 --nsteps 2 -x 0. -S 101 -s 1 -o consume.trees --mutrate 0. --threads
python compare.py
