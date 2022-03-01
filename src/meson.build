
frust_sources = files([
  'main.cpp',
  'frust.cpp',
  'basis.cpp',
  'vertex_data.cpp',
  'write_lattice.cpp',
  'models/cluster_magnet.cpp',
  'models/cluster_magnet_def.cpp',
  'models/model_def.cpp',
  'models/util/lattice.cpp',
])

frust = executable('frust',
  frust_sources,
  dependencies : frust_deps,
  install : true
)
