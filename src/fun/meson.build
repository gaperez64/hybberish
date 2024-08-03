# Library: Functions and expressions
fun_lib = library('fun', files(
    'funexp.c',
    'transformations.c'
    ),
    # IMPORTANT: math functions (floor, ceil, ...)
    # may require explicit linkage to the C math
    # library via the '-lm' gcc flag
    link_args : ['-lm'],
)
