#! /bin/bash
# run:
# bash scripts/build-test-lxplus.sh [noecho] [cleanall] [notest]

# env settings
export LC_CTYPE="C"
export PATH="/afs/cern.ch/user/m/mad/madx/madX:$PATH"

# store lxplus node name in a file
uname -n > build-test-lxplus.run

# error handler
check_error ()
{
	if [ "$?" != "0" ] ; then
		echo -e "\nERROR: $1"
		[ "$2" != "no-exit" ] && exit 1
	fi
}

# I/O redirection
rm -f build-test-lxplus.out
if [ "$1" = "noecho" ] ; then
	shift
	exec > build-test-lxplus.out 2>&1
	check_error "redirection with noecho failed"
else
	exec > >(tee build-test-lxplus.out) 2>&1
	check_error "redirection with tee failed"
fi

echo -e "\n===== Start of build and tests ====="
echo "Date  : `date`"
echo "System: `uname -m -n -r -s`"
echo "Script: $0 $@"

echo -e "\n===== SVN update ====="
svn update
if [ "$?" != "0" ] ; then
	echo -e "\n===== SVN cleanup & update ====="
	svn cleanup
	svn update
	check_error "svn update failed"
fi
# ensure that scripts are executable after an update
chmod u+x scripts/build-test-report.sh $0

echo -e "\n===== Release number ====="
cat VERSION

echo -e "\n===== Clean build ====="
if [ "$1" = "cleanall" ] ; then
	shift
	make cleanall
	check_error "make cleanall failed"
else
	echo "Skipped (no explicit request)."
fi 

echo -e "\n===== Gnu build ====="
source /afs/cern.ch/sw/lcg/contrib/gcc/max/i686-slc6/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux32-gnu
check_error "make all-linux32-gnu failed"

source /afs/cern.ch/sw/lcg/contrib/gcc/max/x86_64-slc6/setup.sh
gcc      --version
g++      --version
gfortran --version
make all-linux64-gnu
check_error "make all-linux64-gnu failed"

echo -e "\n===== Intel build ====="
source /afs/cern.ch/sw/lcg/contrib/gcc/max/i686-slc6/setup.sh
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh ia32
icc      --version
ifort    --version
make all-linux32-intel all-linux32
check_error "make all-linux32-intel failed"

source /afs/cern.ch/sw/lcg/contrib/gcc/max/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/IntelSoftware/linux/all-setup.sh intel64
icc      --version
ifort    --version
make all-linux64-intel all-linux64
check_error "make all-linux64-intel failed"

# license not yet back...
#echo -e "\n===== NagFor build ====="
#export PATH="${PATH}:/afs/cern.ch/sw/fortran/nag2012/bin"
#export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/afs/cern.ch/sw/fortran/nag2012/lib"
#export NAG_KUSARI_FILE="/afs/cern.ch/sw/fortran/nag2012/lib/nag.licence.5.2,lxlic04.cern.ch:"
#make madx-linux-nagfor
#check_error "make madx-linux-nagfor failed" "no-exit"

echo -e "\n===== Lahey 32 build ====="
source /afs/cern.ch/sw/lcg/contrib/gcc/max/i686-slc6/setup.sh
export PATH="${PATH}:/afs/cern.ch/sw/fortran/lahey/lf9562e/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/afs/cern.ch/sw/fortran/lahey/lf9562e/lib"
make madx-linux32-lahey
check_error "make madx-linux32-lahey failed" "no-exit"

echo -e "\n===== Binaries dependencies ====="
make infobindep
check_error "make infobindep failed"

echo -e "\n===== Tests pointless files ====="
make cleantest && make infotestdep
check_error "make infotestdep failed"

echo -e "\n===== Running tests (long) ====="
if [ "$1" = "notest" ] ; then
	shift
	echo "Skipped (explicit request)."
else
	echo ""

	echo -e "\n===== Testing madx-linux64-intel ====="
	make madx-linux64-intel && ls -l madx64 && make cleantest && make tests-all COMP=intel ARCH=64 NOCOLOR=yes
	check_error "make tests-all for madx-linux64-intel failed"

	echo -e "\n===== Testing madx-linux32-intel ====="
	make madx-linux32-intel && ls -l madx32 && make cleantest && make tests-all COMP=intel ARCH=32 NOCOLOR=yes
	check_error "make tests-all for madx-linux32-intel failed"

	echo -e "\n===== Testing madx-linux64-gnu ====="
	make madx-linux64-gnu && ls -l madx64 && make cleantest && make tests-all COMP=gnu ARCH=64 NOCOLOR=yes
	check_error "make tests-all for madx-linux64-gnu failed"

	echo -e "\n===== Testing madx-linux32-gnu ====="
	make madx-linux32-gnu && ls -l madx32 && make cleantest && make tests-all COMP=gnu ARCH=32 NOCOLOR=yes
	check_error "make tests-all for madx-linux32-gnu failed"
fi

# restore the default version
make madx-linux32 > /dev/null && make madx-linux64 > /dev/null
check_error "unable to restore the default version"

# date & end marker
echo -e "\nFinish: `date`"
echo -e "\n===== End of build and tests ====="

# cleanup lxplus node
rm -f build-test-lxplus.run
