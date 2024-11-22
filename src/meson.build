
frust_sources = files([
  'main.cpp',
  'frust.cpp',
  'vertex_data.cpp',
  'write_lattice.cpp',
])

subdir('models')

frust = executable('frust',
  frust_sources,
  dependencies : frust_deps,
  install : true
)
