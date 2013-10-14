# run:
# sh build-test-macosx.sh 1>| build-test-macosx.out 2>&1

echo "===== SVN update ====="
svn update
cat VERSION

echo "===== Gnu build ====="
gcc      --version
g++      --version
gfortran --version
make all-macosx-gnu

echo "===== Intel build ====="
icc      --version
ifort    --version
make all-macosx-intel

echo "===== Dependencies ====="
make infodep

echo "===== Gnu tests (32 bit) ====="
make madx-macosx32-gnu ; make tests-all ARCH=32

echo "===== Gnu tests (64 bit) ====="
make madx-macosx64-gnu ; make tests-all ARCH=64

echo "===== Intel tests (32 bit) ====="
make madx-macosx32-intel ; make tests-all ARCH=32

echo "===== Intel tests (64 bit) ====="
make madx-macosx64-intel ; make tests-all ARCH=64
