project('frust', 'cpp',
  default_options: [
    'buildtype=release',
    'default_library=static',
    'warning_level=3',
    'cpp_std=c++17',
  ]
)

loadleveller_dep = dependency('loadleveller', fallback : [ 'loadleveller', 'loadleveller_dep' ])
eigen_dep = dependency('eigen3', fallback : ['eigen', 'eigen_dep'])
highs_dep = dependency('highs', fallback : [ 'highs', 'highs_dep' ])

downfolded_peierls_coupling = subproject('downfolded_peierls_coupling')
downfolded_peierls_coupling_dep = downfolded_peierls_coupling.get_variable('downfolded_peierls_coupling_dep')

frust_deps = [ loadleveller_dep, eigen_dep, highs_dep, downfolded_peierls_coupling_dep ]

subdir('src')
subdir('test')

