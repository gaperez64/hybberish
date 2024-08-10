# Library: Taylor models and Taylor model flowpipes
taylormodel_lib = library('taylormodel', files(
                            'interval.c',
                            'taylormodel.c',
                            'tmflowpipe.c',
                          ),
                          link_with : [fun_lib, sysode_lib],
                          include_directories : [fun_inc, sysode_inc],
                          # IMPORTANT: math functions (floor, ceil, ...)
                          # may require explicit linkage to the C math
                          # library via the '-lm' gcc flag
                          link_args : ['-lm'],
)
