include(setupNUMDIFFTESTS)
# First parameter is test name,
# second is a bool saying if it is 
# a long test (test-user)
numdiff_test(test-jacobian False)
numdiff_test(test-jacobian-2 False)
numdiff_test(test-jacobian-knobs False)
numdiff_test(test-ptc-twiss False)
