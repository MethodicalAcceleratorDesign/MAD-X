# run:
# sh scripts/build-test-lxplus.sh 1> build-test-lxplus.out 2>&1
# tail -f build-test-lxplus.out

export LC_CTYPE="C"

echo -e "\n===== Start of build and tests ====="
date
uname -m -n -r -s

echo -e "\n===== SVN update ====="
svn update && cat VERSION

echo -e "\n===== Clean build ====="
make cleanall ; make cleanall ARCH=32

echo -e "\n===== Gnu build ====="
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/i686-slc6-gcc48-opt/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux32-gnu
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/setup.sh
make all-linux64-gnu

echo -e "\n===== Intel build ====="
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh ia32
icc      --version
ifort    --version
make all-linux32-intel all-linux32
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh intel64
make all-linux64-intel all-linux64

echo -e "\n===== Dependencies ====="
make infodep

echo -e "\n===== Gnu tests (32 bit) ====="
make madx-linux32-gnu && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes

echo -e "\n===== Gnu tests (64 bit) ====="
make madx-linux64-gnu && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes

echo -e "\n===== Intel tests (32 bit) ====="
make madx-linux32-intel && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes

echo -e "\n===== Intel tests (64 bit) ====="
make madx-linux64-intel && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes

echo -e "\n===== End of build and tests ====="
date
