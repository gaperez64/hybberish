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
               link_with : taylormodel_lib,
               include_directories : taylormodel_inc,
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
