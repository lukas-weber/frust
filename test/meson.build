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

sim_tests = find_program('simtest.py')
test('ed_square',
   sim_tests,
   args : [frust, files('testjobs/square')],
   timeout : 200,
   is_parallel : false,
   suite : 'ed_test')

test('ed_triangle',
   sim_tests,
   args : [frust, files('testjobs/triangle')],
   timeout : 200,
   is_parallel : false,
   suite : 'ed_test')

test('ed_bilayer',
   sim_tests,
   args : [frust, files('testjobs/bilayer')],
   timeout : 200,
   is_parallel : false,
   suite : 'ed_test')

test('ed_dimerized_bilayer',
   sim_tests,
   args : [frust, files('testjobs/dimerized_bilayer')],
   timeout : 200,
   is_parallel : false,
   suite : 'ed_test')

test('ed_shastry_sutherland',
   sim_tests,
   args : [frust, files('testjobs/shastry_sutherland')],
   timeout : 200,
   is_parallel : false,
   suite : 'ed_test')

test('ed_triang_square',
   sim_tests,
   args : [frust, files('testjobs/triangle_square')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

test('ed_triang_square_rot',
   sim_tests,
   args : [frust, files('testjobs/triangle_square_rot')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

test('ed_triang_square3',
   sim_tests,
   args : [frust, files('testjobs/triangle_square3')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

test('ed_diamond_square_dim_spi',
   sim_tests,
   args : [frust, files('testjobs/diamond_square_dim_spi')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

# FIXME: this is a test for the rotated dimer basis. currently JDim does not
# compare correctly to the ed results which do not know about the basis rotation
#
#test('ed_triang_square4',
#   sim_tests,
#   args : [frust, files('testjobs/triangle_square4')],
#   timeout : 5000,
#   is_parallel : false,
#   suite : 'ed_test')

# XXX: maxwormlen is not confirmed to work. I am leaving the functionality enabled for now
# in case somebody needs it. Still, I am not sure if it works correctly.
#
#test('ed_maxwormlen',
#   sim_tests,
#   args : [frust, files('testjobs/maxwormlen')],
#   timeout : 5000,
#   is_parallel : false,
#   suite : 'ed_test')

test('ed_chirality',
   sim_tests,
   args : [frust, files('testjobs/chirality')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

test('ed_kagome',
   sim_tests,
   args : [frust, files('testjobs/kagome')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

test('ed_kagome_tri',
   sim_tests,
   args : [frust, files('testjobs/kagome_tri')],
   timeout : 5000,
   is_parallel : false,
   suite : 'ed_test')

#test('ed_kagome_dim_spi',
#   sim_tests,
#   args : [frust, files('testjobs/kagome_dim_spi')],
#   timeout : 5000,
#   is_parallel : false,
#   suite : 'ed_test')

hamiltonian_test = find_program('test_equal_hamiltonian.py')
test('equal_hamiltonian', hamiltonian_test, args : [frust, files('testjobs/triangle_square_bases')])
