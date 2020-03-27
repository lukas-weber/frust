project('frust', 'cpp',
  default_options: [
    'buildtype=release',
    'default_library=static',
    'warning_level=3',
    'cpp_std=c++17',
  ]
)

loadleveller_dep = dependency('loadleveller', fallback : [ 'loadleveller', 'loadleveller_dep' ])
solp_dep = dependency('solp', fallback : [ 'solp', 'solp_dep' ])
catch2_dep = dependency('catch2', fallback : [ 'catch2', 'catch2_dep' ])

frust_deps = [ loadleveller_dep, solp_dep ]

subdir('src')
#subdir('test')
