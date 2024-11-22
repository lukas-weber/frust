test_srcs = files([
  'test_opercode.cpp',
#  'test_wormfunc.cpp'
])

catch2_dep = dependency('catch2', fallback : [ 'catch2', 'catch2_dep' ])
unit_tests = executable('unit_tests', test_srcs,
                        dependencies : [ frust_deps, catch2_dep ],
			include_directories : include_directories('../src'))

test('unit tests', unit_tests)
