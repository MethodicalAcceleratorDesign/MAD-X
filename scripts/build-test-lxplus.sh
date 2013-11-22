# run:
# sh scripts/build-test-lxplus.sh [cleanall]
# tail -f build-test-lxplus.out

# env settings
export LC_CTYPE="C"
export PATH="/afs/cern.ch/user/m/mad/madx/madX:$PATH"

# I/O redirection
rm -f build-test-lxplus.out
exec 1> build-test-lxplus.out 2>&1
uname -n > build-test-lxplus.run

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
	make cleanall && make cleanall ARCH=32
	[ "$?" != "0" ] && echo "ERROR: make cleanall failed"
else
	echo "Skipped (no explicit request)."
fi 

echo -e "\n===== Gnu build ====="
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/i686-slc6-gcc48-opt/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux32-gnu
# to handle bad fortran compiler, restart from scratch (only once)
if [ "$?" != "0" ] ; then 
	make cleanall && make cleanall ARCH=32
	[ "$?" != "0" ] && echo "ERROR: make cleanall failed"
	make all-linux32-gnu
	[ "$?" != "0" ] && echo "ERROR: make all-linux32-gnu failed"
fi

source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux64-gnu
[ "$?" != "0" ] && echo "ERROR: make all-linux64-gnu failed"

echo -e "\n===== Intel build ====="
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh ia32
icc      --version
ifort    --version
make all-linux32-intel all-linux32
[ "$?" != "0" ] && echo "ERROR: make all-linux32-intel failed"

source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh intel64
icc      --version
ifort    --version
make all-linux64-intel all-linux64
[ "$?" != "0" ] && echo "ERROR: make all-linux64-intel failed"

echo -e "\n===== Binaries dependencies ====="
make infobindep
[ "$?" != "0" ] && echo "ERROR: make infobindep failed"

echo -e "\n===== Tests pointless files ====="
make cleantest && make infotestdep
[ "$?" != "0" ] && echo "ERROR: make infotestdep failed"

echo -e "\n===== Testing madx-linux64-intel ====="
make madx-linux64-intel && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux64-intel failed"

echo -e "\n===== Testing madx-linux32-intel ====="
make madx-linux32-intel && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux32-intel failed"

echo -e "\n===== Testing madx-linux64-gnu ====="
make madx-linux64-gnu && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux64-gnu failed"

echo -e "\n===== Testing madx-linux32-gnu ====="
make madx-linux32-gnu && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-linux32-gnu failed"

# restore the default version
make madx-linux32 > /dev/null && make madx-linux64 > /dev/null
[ "$?" != "0" ] && echo "ERROR: error restoring the default version"

# cleanup lxplus node
rm -f build-test-lxplus.run

# date & end marker
date
echo -e "\n===== End of build and tests ====="
