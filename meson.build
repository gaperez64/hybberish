project('hybberish', 'c',
  version : '0.1',
  default_options : ['c_std=c11',
                     'warning_level=3']
)

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
                 arguments: ['-Wall',
                             '--output=@OUTPUT0@',
                             '--defines=@OUTPUT1@',
                             '@INPUT@'])

# Tests

#exe = executable('hybberish', 'hybberish.c',
#  install : true)

# test('basic', exe)
