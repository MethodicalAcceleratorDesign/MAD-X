# run:
# sh scripts/build-test-macosx.sh 1> build-test-macosx.out 2>&1
# tail -f build-test-macosx.out

export LC_CTYPE="C"
export PATH=/opt/local/bin:$PATH

echo "\n===== Start of build and tests ====="
date
uname -n -r -s -m

echo "\n===== SVN update ====="
svn update && cat VERSION

echo "\n===== Clean build ====="
make cleanall ; make cleanall ARCH=32

echo "\n===== Gnu build ====="
gcc      --version
g++      --version
gfortran --version
make all-macosx-gnu

echo "\n===== Intel build ====="
icc      --version
ifort    --version
make all-macosx-intel all-macosx

echo "\n===== Dependencies ====="
make infodep

echo "\n===== Gnu tests (32 bit) ====="
make madx-macosx32-gnu && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes

echo "\n===== Gnu tests (64 bit) ====="
make madx-macosx64-gnu && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes

echo "\n===== Intel tests (32 bit) ====="
make madx-macosx32-intel && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes

echo "\n===== Intel tests (64 bit) ====="
make madx-macosx64-intel && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes

echo "\n===== End of build and tests ====="
date
