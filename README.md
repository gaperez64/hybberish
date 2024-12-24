# Hybberish
Utilities for the anylisis of hybrid systems

## Notes on how to compile

The compilation is quite straightforward.

A note for MacOS users: the "bison" utility that ships with MacOS is not compatible. Please install a newer version using homebrew. GNU Bison 3.8.2 should work.

## Notes on generating documentation

[Doxygen](https://www.doxygen.nl/manual/index.html) is used for documentation.

To generate the documentation locally, follow the following steps (on linux), starting at root of the project.

```sh
# Have doxygen installed
apt install doxygen

# Generate documentation in the doxygen dir
cd doxygen/
doxygen
```

## Running the notebooks

There are a number of Pluto notebooks present in this repository. Pluto can be started using `./run_pluto.sh`.

The dependency situation is a bit complicated. This project requires the latest version of the TaylorSeries.jl package (v0.18.2), but the latest version of TaylorModels.jl only supports versions 0.17.x.

The dependencies can be resolved by opening Pkg-mode and activating the current environment. The commands `instantiate` and `rm TaylorModels` should be ran. Next, clone the reposity of the TaylorModels.jl package in a separate directory. Modify its `Project.toml` file and set the required version of TaylorSeries.jl to `=0.18` instead of `=0.17`. Then, in the Pkg-mode of the Hybberish project, develop your local modification of TaylorModels.jl by issuing the command `develop <PATH TO LOCAL TAYLOR MODELS JL PACKGE>`. The notebooks should now run normally.



