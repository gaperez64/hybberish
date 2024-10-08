# Flex as a generator
flex = find_program('flex')
lgen = generator(flex,
                 output: ['@BASENAME@.lex.c', '@BASENAME@.lex.h'],
                 arguments: ['--outfile=@OUTPUT0@',
                             '--header-file=@OUTPUT1@',
                             '@INPUT@'])

# Bison as a generator
bison = find_program('bison')
pgen = generator(bison,
                 output: ['@BASENAME@.tab.c', '@BASENAME@.tab.h'],
                 arguments: ['--output=@OUTPUT0@',
                             '--defines=@OUTPUT1@',
                             '@INPUT@'])

# ODE parsing library
odeparse_lib = library('odeparse',
                       [lgen.process('odes.l'),
                        pgen.process('odes.y')],
                       link_with : [fun_lib, sysode_lib],
                       include_directories : [fun_inc, sysode_inc])

# variable valuation parsing library
varparse_lib = library('varparse',
                       [lgen.process('vars.l'),
                        pgen.process('vars.y')],
                       link_with : [fun_lib, varmath_lib],
                       include_directories : [fun_inc, varmath_inc])
