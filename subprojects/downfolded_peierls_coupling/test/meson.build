test_srcs = files(
  'unit_tests.cpp'
)

catch2_dep = dependency('catch2', fallback : [ 'catch2', 'catch2_dep' ])
unit_tests = executable('unit_tests', test_srcs,
                        dependencies : [ downfolded_peierls_coupling_dep, catch2_dep ],
			include_directories : include_directories('../src'))

test('unit_tests', unit_tests)
