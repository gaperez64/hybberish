project('hybberish', 'c',
  version : '0.1',
  default_options : ['warning_level=3']
)

subdir('src/fun')
fun_inc = include_directories('src/fun')
subdir('src/systems')
sysode_inc = include_directories('src/systems')
subdir('src/parsers')
odeparse_inc = include_directories('src/parsers')

subdir('tests')
