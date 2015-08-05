#!/usr/bin/env python

from __future__ import print_function


import sys,os,subprocess


# Find tests defined in makefile:

fin=open('Makefile_test')
lines=[]
nextline=''
for l in fin:
    if len(l.strip())>0 and l.strip()[-1]=='\\':
        nextline+=l.strip()[:-1]+' '
    else:
        nextline+=l.strip()
        if len(nextline)>0:
            lines.append(nextline)
        nextline=''

for l in lines:
    if l.split(':=')[0].strip()=='tests-all':
        tests_all=l.split()[2:]
    if l.split(':=')[0].strip()=='tests-long':
        tests_long=l.split()[2:]

# Find tests defined in cmake:

if not os.path.isfile('build/CMakeCache.txt'):
    raise ValueError('Need a build directory configured..')

os.chdir('build')

# Python3 equivalent:
# ctest_stdout=str(subprocess.check_output(['ctest','-N']),encoding='utf8').split('\n')
ctest_stdout=subprocess.check_output(['ctest','-N']).split('\n')
ctest_all=[]
ctest_long=[]
for l in ctest_stdout:
    lsp=l.split()
    if len(lsp)==3 and lsp[0]=='Test' and '#' in lsp[1]:
        if lsp[2][-4:]=='LONG':
            test_name=lsp[2][:-5]
            ctest_long.append(test_name)
        else:
            test_name=lsp[2]
        ctest_all.append(test_name)


# Compare the lists:

ctest_missing=[]
ctest_missing_long=[]
long_only_make=[]
long_only_cmake=[]
for test in tests_all:
    if test[:4]!='test': # ignore difficult make syntax for now..
        continue
    if test not in ctest_all:
        ctest_missing.append(test)
    if test in tests_long and not test in ctest_long:
        ctest_missing_long.append(test)
        long_only_make.append(test)
    if test in ctest_long and not test in tests_long:
        long_only_cmake.append(test)

if long_only_make:
    print("Tests that are only long in make:")
    print(" ".join(long_only_make))
if long_only_cmake:
    print("Tests that are only long in cmake:")
    print(" ".join(long_only_cmake))

# Print the needed additions:
if ctest_missing:
    print(" -- Missing tests: --")
    for test in ctest_missing:
        if test in ctest_missing_long:
            print('numdiff_test('+test+' 1)')
        else:
            print('numdiff_test('+test+' 0)')

# Print tests that don't exist:
for test in ctest_all:
    path=os.path.join('..','tests',test,test+'.madx')
    if not os.path.exists(path):
        print("{} is missing".format(test))
