# run:
# sh build-test-macosx.sh 1>| build-test-macosx.out 2>&1
# tail -f build-test-macosx.out

echo "\n===== SVN update ====="
svn update
cat VERSION

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
make madx-macosx32-gnu && ls -l madx32 && make tests-all ARCH=32

echo "\n===== Gnu tests (64 bit) ====="
make madx-macosx64-gnu && ls -l madx64 && make tests-all ARCH=64

echo "\n===== Intel tests (32 bit) ====="
make madx-macosx32-intel && ls -l madx32 && make tests-all ARCH=32

echo "\n===== Intel tests (64 bit) ====="
make madx-macosx64-intel && ls -l madx64 && make tests-all ARCH=64
