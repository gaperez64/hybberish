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
