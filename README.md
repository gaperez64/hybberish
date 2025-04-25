# Hybberish
Utilities for the analysis of hybrid systems

## Dependencies

The Julia dependencies are listed in the root `project.toml` file. The **specific versions** of *TaylorSeries* and *IntervalArithmetic* **listed there must be installed**, instead of the most recent versions. We started development making use of the *TaylorModels* package, which used an older version of *IntervalArithmetic*. The use of *IntervalBox* now locks us into an older version of *IntervalArithmetic*, which consequently limits the version of *TaylorSeries* as well, until we refactor the code to allow bumping the *IntervalArithmetic* version.

## How to run

The project is currently unpolished and untested. Thus, the examples contained in `examples/` are the most interesting way to gauge the performance of the current implementation.

Once all dependencies are installed, simply invoking the example scripts should suffice to generate a plot of the results. The scripts that are currently known to at least run, and which are interesting, are listed below. Note that some of the examples may take a few minutes to finish running.

```sh
# Invoke an example script.
julia examples/simple_trucktrailer.jl
julia examples/example_3_3_10.jl
julia examples/example_m_neher.jl
julia examples/simple_example.jl
julia examples/simple_example_modified.jl
julia examples/example_tmjets.jl
```


## Notes on Implementation

Confusingly, the file `src/tm_integration_qr.jl` implements our most up-to-date version of naive Taylor model integration. The original intent was to add QR preconditioning to improve the naive method, but this was not finished on time.

The file `src/tm_integration.jl` implements an implementation of naive Taylor model integration that we want to deprecate and replace by the newer version. Nevertheless, this older version was used to pad the NN trace for the DAggr approach.
