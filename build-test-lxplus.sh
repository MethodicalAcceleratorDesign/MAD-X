# run:
# sh build-test-lxplus.sh 1>| build-test-lxplus.out 2>&1

echo "===== SVN update ====="
svn update
cat VERSION

echo "===== Gnu build ====="
gcc      --version
g++      --version
gfortran --version
make all-linux-gnu 

echo "===== Intel build ====="
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh ia32
icc      --version
ifort    --version
make all-linux32-intel
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh intel64
make all-linux64-intel

echo "===== Dependencies ====="
make infodep

echo "===== Gnu tests (32 bit) ====="
make madx-linux32-gnu ; make tests-all ARCH=32

echo "===== Gnu tests (64 bit) ====="
make madx-linux64-gnu ; make tests-all ARCH=64

echo "===== Intel tests (32 bit) ====="
make madx-linux32-intel ; make tests-all ARCH=32

echo "===== Intel tests (64 bit) ====="
make madx-linux64-intel ; make tests-all ARCH=64
