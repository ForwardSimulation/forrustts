#!/bin/bash

# ./target/release/examples/forward_simulation --popsize 5000 --nsteps 10000 -x 1e-3 -S 101 -s 100 -o regular.trees --mutrate 0.
# ./target/release/examples/forward_simulation --popsize 5000 --nsteps 10000 -x 1e-3 -S 101 -s 100 -o consume.trees --mutrate 0. --threads
./target/release/examples/forward_simulation --popsize 2 --nsteps 2 -x 0. -S 101 -s 1 -o regular.trees --mutrate 0.
./target/release/examples/forward_simulation --popsize 2 --nsteps 2 -x 0. -S 101 -s 1 -o consume.trees --mutrate 0. --threads
python compare.py
