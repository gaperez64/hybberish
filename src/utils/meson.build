# Library: Utility functions
utils_lib = library('utils', files(
                      'utils.c',
                    ),
                    # IMPORTANT: math functions (floor, ceil, ...)
                    # may require explicit linkage to the C math
                    # library via the '-lm' gcc flag
                    link_args : ['-lm'],
)
