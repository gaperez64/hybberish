project('hybberish', 'c',
  version : '0.1',
  default_options : ['c_std=c23',
                     'warning_level=3']
)

subdir('src/parsers')
subdir('src/systems')
sysode_inc = include_directories('src/systems')

subdir('tests')
