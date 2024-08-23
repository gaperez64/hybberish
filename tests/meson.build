# Tests
t = executable('fun_exp_test', 'fun_exp_test.c',
               link_with : fun_lib,
               include_directories : fun_inc)
test('test sysode expressions', t)

t = executable('sysode_explist_test', 'sysode_explist_test.c',
               link_with : [fun_lib, sysode_lib],
               include_directories : [fun_inc, sysode_inc])
test('test sysode expression lists', t)

t = executable('odeparse_test', 'odeparse_test.c',
               link_with : [fun_lib, sysode_lib, odeparse_lib],
               include_directories : [fun_inc, sysode_inc, odeparse_inc])
test('test ode parser', t)

t = executable('derivative_test', 'derivative_test.c',
               link_with : fun_lib,
               include_directories : fun_inc)
test('test derivative', t)

t = executable('integral_test', 'integral_test.c',
               link_with : fun_lib,
               include_directories : fun_inc)
test('test integral', t)

t = executable('interval_test', 'interval_test.c',
               link_with : varmath_lib,
               include_directories : varmath_inc,
               # IMPORTANT: math functions (floor, ceil, ...)
               # may require explicit linkage to the C math
               # library via the '-lm' gcc flag
               link_args : ['-lm'],
               )
test('test interval', t)


t = executable('transformations_test', 'transformations_test.c',
               link_with : fun_lib,
               include_directories : fun_inc,
               )
test('test expression tree transformations', t)

t = executable('fun_exp_equivalence_test', 'fun_exp_equivalence_test.c',
               link_with : fun_lib,
               include_directories : fun_inc,
               )
test('test expression tree equivalence functions', t)

t = executable('fun_exp_polynomials_test', 'fun_exp_polynomials_test.c',
               link_with : fun_lib,
               include_directories : fun_inc,
               # IMPORTANT: math functions (floor, ceil, ...)
               # may require explicit linkage to the C math
               # library via the '-lm' gcc flag
               link_args : ['-lm'],
               )
test('test fun exp functions that operate on polynomials', t)

t = executable('definite_integral_test', 'definite_integral_test.c',
               link_with : fun_lib,
               include_directories : fun_inc)
test('test definite integral', t)

t = executable('tmflowpipe_test', 'tmflowpipe_test.c',
               link_with : [fun_lib, varmath_lib, sysode_lib, taylormodel_lib],
               include_directories : [fun_inc, varmath_inc, sysode_inc, taylormodel_inc],
               # IMPORTANT: math functions (floor, ceil, ...)
               # may require explicit linkage to the C math
               # library via the '-lm' gcc flag
               link_args : ['-lm'],
               )
test('test taylor model flowpipe overapprox computations', t)

t = executable('taylormodel_test', 'taylormodel_test.c',
               link_with : [fun_lib, varmath_lib, sysode_lib, taylormodel_lib],
               include_directories : [fun_inc, varmath_inc, taylormodel_inc],
               )
test('test taylor model interface', t)

t = executable('varparse_test', 'varparse_test.c',
               link_with : [fun_lib, varmath_lib, varparse_lib],
               include_directories : [fun_inc, varmath_inc, odeparse_inc])
test('test var parser', t)
