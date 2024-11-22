project('downfolded_peierls_coupling', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++17'])

subdir('src')

libdownfolded_peierls_coupling = library('downfolded_peierls_coupling',
  downfolded_peierls_coupling_sources)

install_headers(downfolded_peierls_coupling_headers)

downfolded_peierls_coupling_dep = declare_dependency(
  include_directories : include_directories('./src'),
  link_with : libdownfolded_peierls_coupling,
  version : meson.project_version())

pymod = import('python')
py = pymod.find_installation()

downfolded_peierls_coupling_py = files('__init__.py') 


py.extension_module('downfolded_peierls_coupling',
  downfolded_peierls_coupling_sources,
  install: true,
  subdir : 'downfolded_peierls_coupling'
)

py.install_sources(downfolded_peierls_coupling_py + downfolded_peierls_coupling_headers,
  pure : false,
  subdir : 'downfolded_peierls_coupling'
)

if not meson.is_subproject()
  subdir('test')
endif
