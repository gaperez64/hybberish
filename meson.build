project('hybberish', 'c',
  version : '0.1',
  default_options : ['warning_level=3']
)

subdir('src/fun')
fun_inc = include_directories('src/fun')
subdir('src/systems')
sysode_inc = include_directories('src/systems')
subdir('src/varmath')
varmath_inc = include_directories('src/varmath')
subdir('src/parsers')
odeparse_inc = include_directories('src/parsers')
subdir('src/taylormodel')
taylormodel_inc = include_directories('src/taylormodel')

subdir('tests')
