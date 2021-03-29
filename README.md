# Forward simulation with tree sequence recording in rust

This package is currently "experimental"!

`forrustts` (pronounced "forests") is a port of many ideas from the [fwdpp](https://github.com/molpopgen/fwdpp) library from C++ to rust.

It is licensed under the GNU General Public License, version 3 or later ("GPL3+").

## Getting started

Install the GNU Scientific Library.
For example:

```
apt install libgsl-dev
```

You may use `conda` or `brew` as you see fit.

[Install rust](https://www.rust-lang.org/learn/get-started)

Then:

```
cargo build
cargo test
```

Example programs are in the subdirectory `forrustts_examples`:

```
cd forrustts_examples
cargo build --release
```

The binaries will then be found in `target/release/`.

These programs use [clap](https://crates.io/crates/clap) for command-line options.
Pass ``--help`` to any of them for usage information.

## Development information

## CI

CI testing is done using GitHub actions for both `Linux` and `macOS`.
These actions include using [clippy](https://crates.io/crates/clippy/0.0.211), which is a very strict code linter.
The actions also check code format using [rustfmt](https://crates.io/crates/rustfmt-nightly).

### Code coverage

Use [tarpaulin](https://docs.rs/crate/cargo-tarpaulin/0.3.12).
The documentation for that crate is excellent.
The short version is:

```
cargo tarpaulin -o html
```

This command will run the tests and generate a nice `html` report.

## Change log

## 2021-03-29, Version 0.1.3

#### Commits
- [[`99821a933c`](https://github.com/molpopgen/forrustts/commit/99821a933c473f77e483b4929992dfeb9e445680)] Bump version to 0.1.3 (molpopgen)
- [[`a34d582e02`](https://github.com/molpopgen/forrustts/commit/a34d582e027a0992939ae41206d3c87204d30a0b)] add change log info for 0.1.2 (molpopgen)
- [[`cfcfb7e657`](https://github.com/molpopgen/forrustts/commit/cfcfb7e6574a50abd12d77cd71aa824d224b4f23)] Merge pull request #72 from molpopgen/tskit_tools (Kevin R. Thornton)
- [[`d69bf0889e`](https://github.com/molpopgen/forrustts/commit/d69bf0889eb8b5e62a4c36c14ac384ec3524d812)] update name in examples (molpopgen)
- [[`7e65155553`](https://github.com/molpopgen/forrustts/commit/7e651555532e0a26373ac61e2f5c329fcb577f0a)] rename crate::tskit to crate::tskit_tools (molpopgen)
- [[`dffc91a41c`](https://github.com/molpopgen/forrustts/commit/dffc91a41ca23a244500c324c615da022b3460c7)] Merge pull request #71 from molpopgen/fix_doc_issues (Kevin R. Thornton)
- [[`a6b929818c`](https://github.com/molpopgen/forrustts/commit/a6b929818cafcb88f87344eba7aba3ec6a382967)] Fix various links in docs. (molpopgen)
- [[`9e56d40662`](https://github.com/molpopgen/forrustts/commit/9e56d40662c900fc50a5ed422d3f87a6ced106a5)] Merge pull request #69 from molpopgen/streamline_simplification_code_organization (Kevin R. Thornton)
- [[`e441f9202a`](https://github.com/molpopgen/forrustts/commit/e441f9202a49f41ad37613461df8c970ee06c131)] Move all simplification code into one source file (molpopgen)
- [[`4123c16adf`](https://github.com/molpopgen/forrustts/commit/4123c16adfc71f0763a8545d3b6f8e57fd4ae3f0)] Merge pull request #70 from molpopgen/update_cancel_action (Kevin R. Thornton)
- [[`f694d05abd`](https://github.com/molpopgen/forrustts/commit/f694d05abdbe9b745d83d534dd97b85f078c922d)] Add cancel to macos.  Update token. (molpopgen)
- [[`22150f32fe`](https://github.com/molpopgen/forrustts/commit/22150f32fe242f0042962eadfb4992ff81a40abb)] Merge pull request #68 from ForwardSimulation/dependabot/add-v2-config-file (Kevin R. Thornton)
- [[`3c9ed4f211`](https://github.com/molpopgen/forrustts/commit/3c9ed4f21127c357406cce31bac5e757f71adfef)] Create Dependabot config file (dependabot-preview[bot])
- [[`01a529220c`](https://github.com/molpopgen/forrustts/commit/01a529220cd5e4e1523d8bb19d5ca7a71d510b5d)] Merge pull request #67 from ForwardSimulation/dependabot/cargo/GSL-4.0.0 (dependabot-preview[bot])
- [[`c8567deedf`](https://github.com/molpopgen/forrustts/commit/c8567deedfc3c64babd6d8e57596e6ceb88fbddf)] Update GSL requirement from 2.0.1 to 4.0.0 (dependabot-preview[bot])


### 2021-03-26, Version 0.1.2

#### Commits
- [[`21a77982a8`](https://github.com/molpopgen/forrustts/commit/21a77982a8530ae32791ea7eee35af8b9f6da0a7)] bump version to 0.1.2 (molpopogen)
- [[`4653a8b5d1`](https://github.com/molpopgen/forrustts/commit/4653a8b5d128405e9da6f50c58caadb92d97beca)] Merge pull request #64 from molpopgen/update_tskit_dependency (Kevin R. Thornton)
- [[`6bc5030d9e`](https://github.com/molpopgen/forrustts/commit/6bc5030d9e349a0d059f1d06cb741a94522cee01)] Fix string format lint in assert. (molpopogen)
- [[`f36c6c91df`](https://github.com/molpopgen/forrustts/commit/f36c6c91df6500ecc6f2c31246bcf5d5b32bd59e)] update to tskit 0.1.1 (molpopogen)
- [[`b1424be2f8`](https://github.com/molpopgen/forrustts/commit/b1424be2f80400c799503c08500cf48b86c6cfd1)] Merge pull request #62 from molpopgen/doc_fixes (Kevin R. Thornton)
- [[`5771eb18d6`](https://github.com/molpopgen/forrustts/commit/5771eb18d6659af55912697bb7bfeac2b76f2c48)] Document other differences from tskit. (molpopgen)


