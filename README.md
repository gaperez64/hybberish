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
