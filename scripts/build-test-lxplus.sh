# run:
# sh scripts/build-test-lxplus.sh [cleanall]
# tail -f build-test-lxplus.out

# I/O redirection
rm -f build-test-lxplus.out
exec 1> build-test-lxplus.out 2>&1

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

echo -e "\n===== Start of build and tests ====="
date
uname -m -n -r -s

echo -e "\n===== SVN update ====="
svn update
[ "$?" != "0" ] && echo "ERROR: svn update failed"

echo -e "\n===== Release number ====="
cat VERSION

echo -e "\n===== Clean build ====="
if [ "$1" = "cleanall" ] ; then
	make cleanall SHOW=yes ; make cleanall ARCH=32 SHOW=yes
else
	echo "not requested, skipped"
fi 

echo -e "\n===== Gnu build ====="
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/i686-slc6-gcc48-opt/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux32-gnu SHOW=yes
[ "$?" != "0" ] && echo "ERROR: make all-linux32-gnu failed"
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/setup.sh
make all-linux64-gnu SHOW=yes
[ "$?" != "0" ] && echo "ERROR: make all-linux64-gnu failed"

echo -e "\n===== Intel build ====="
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh ia32
icc      --version
ifort    --version
make all-linux32-intel all-linux32 SHOW=yes
[ "$?" != "0" ] && echo "ERROR: make all-linux32-intel failed"
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh intel64
make all-linux64-intel all-linux64 SHOW=yes
[ "$?" != "0" ] && echo "ERROR: make all-linux64-intel failed"

echo -e "\n===== Dependencies ====="
make infodep SHOW=yes

echo -e "\n===== Gnu tests (32 bit) ====="
make madx-linux32-gnu SHOW=yes && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux32-gnu failed"

echo -e "\n===== Gnu tests (64 bit) ====="
make madx-linux64-gnu SHOW=yes && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux64-gnu failed"

echo -e "\n===== Intel tests (32 bit) ====="
make madx-linux32-intel SHOW=yes && ls -l madx32 && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux32-intel failed"

echo -e "\n===== Intel tests (64 bit) ====="
make madx-linux64-intel SHOW=yes && ls -l madx64 && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux64-intel failed"

echo -e "\n===== End of build and tests ====="
date
