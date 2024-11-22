test_srcs = files([
  'test_opercode.cpp',
  'test_util.cpp',
  '../src/basis.cpp'
])

catch2_dep = dependency('catch2', fallback : [ 'catch2', 'catch2_dep' ])
unit_tests = executable('unit_tests', test_srcs,
                        dependencies : [ frust_deps, catch2_dep ],
			include_directories : include_directories('../src'))

test('unit_tests', unit_tests)

subdir('testjobs')

fs = import('fs')

test_simulation = find_program('test_simulation.py')
gen_results = find_program('gen_result_files')

foreach tj : test_jobfiles
  test('ed_' + fs.name(tj),
    test_simulation,
    args : [frust, tj, '--result-file', join_paths(meson.project_source_root(), 'testjobs/ed/', fs.name(tj) + '.ed.json')],
    timeout : 5000,
    is_parallel : false,
    suite : 'ed')
  test('seeded_' + fs.name(tj),
    test_simulation,
    timeout : 5000,
    args : [frust, tj, '--seeded', '--result-file', join_paths(meson.project_source_root(), 'test/testjobs/seeded/', fs.name(tj) + '.seeded.json')],
    suite : 'seeded')
endforeach

run_target('gen-seeded-results',
  command : [gen_results, frust, 'seeded', test_jobfiles])

hamiltonian_test = find_program('test_equal_hamiltonian.py')
test('equal_hamiltonian', hamiltonian_test, args : [frust, files('testjobs/triangle_square_bases')])
