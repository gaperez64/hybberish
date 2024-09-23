# Library: Variable valuations and interval math
varmath_lib = library('varmath', files(
                        'interval.c',
                        'variables.c',
                      ),
                      # IMPORTANT: math functions (floor, ceil, ...)
                      # may require explicit linkage to the C math
                      # library via the '-lm' gcc flag
                      link_args : ['-lm'],
)
