# Tests
t = executable('sysode_exp_test', 'sysode_exp_test.c',
               link_with : sysode_lib,
               include_directories : sysode_inc)
test('test sysode expressions', t)

t = executable('sysode_explist_test', 'sysode_explist_test.c',
               link_with : sysode_lib,
               include_directories : sysode_inc)
test('test sysode expression lists', t)

t = executable('odeparse_test', 'odeparse_test.c',
               link_with : [sysode_lib, odeparse_lib],
               include_directories : [sysode_inc, odeparse_inc])
test('test ode parser', t)
