# Library: Taylor models and Taylor model flowpipes
taylormodel_lib = library('taylormodel', 'interval.c',
                          # IMPORTANT: math functions (floor, ceil, ...)
                          # may require explicit linkage to the C math
                          # library via the '-lm' gcc flag
                          link_args : ['-lm'],
)
