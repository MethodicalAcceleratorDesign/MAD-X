#! /bin/bash
# run:
# bash scripts/build-test-linux.sh [noecho] [cleanall] [notest]

# env settings
export PATH="`pwd`:/opt/intel/bin:$PATH"

# error handler
check_error ()
{
	if [ "$?" != "0" ] ; then
		echo -e "\nERROR: $1"
		[ "$2" != "no-exit" ] && exit 1
	fi
}

# I/O redirection
rm -f build-test-linux.out
if [ "$1" = "noecho" ] ; then
	shift
	export NOCOLOR=yes
	exec > build-test-linux.out 2>&1
	check_error "redirection with noecho failed"
else
	exec > >(tee build-test-linux.out) 2>&1
	check_error "redirection with tee failed"
fi

echo -e "\n===== Start of build and tests ====="
echo "Date  : `date`"
echo "UserId: `whoami`"
echo "System: `uname -m -n -r -s`"
echo "Script: $0 $@"
echo "Shell : `echo $SHELL`"
echo "User  : `echo $USER`"
echo "Home  : `echo $HOME`"
echo "Lang  : `echo $LANG`"

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
gcc      --version
g++      --version
gfortran --version

make all-linux32-gnu
check_error "make all-linux32-gnu failed" "no-exit"

make all-linux64-gnu
check_error "make all-linux64-gnu failed" "no-exit"

echo -e "\n===== Intel build ====="

source compilervars.sh ia32
icc      -V
ifort    -V
make all-linux32-intel
check_error "make all-linux32-intel failed" "no-exit"

source compilervars.sh intel64
icc      -V
ifort    -V
make all-linux64-intel
check_error "make all-linux64-intel failed" "no-exit"

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

	echo -e "\n===== Testing madx-linux32-intel ====="
	make madx-linux32-intel && ls -l madx32 && make cleantest && make tests-all COMP=intel ARCH=32 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-linux32-intel failed" "no-exit"

	echo -e "\n===== Testing madx-linux64-intel ====="
	make madx-linux64-intel && ls -l madx64 && make cleantest && make tests-all COMP=intel ARCH=64 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-linux64-intel failed" "no-exit"

	echo -e "\n===== Testing madx-linux32-gnu ====="
	make madx-linux32-gnu && ls -l madx32 && make cleantest && make tests-all COMP=gnu ARCH=32 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-linux32-gnu failed" "no-exit"

	echo -e "\n===== Testing madx-linux64-gnu ====="
	make madx-linux64-gnu && ls -l madx64 && make cleantest && make tests-all COMP=gnu ARCH=64 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-linux64-gnu failed" "no-exit"
fi

# restore the default version
make madx-linux32-gnu > /dev/null && make madx-linux64-gnu > /dev/null
check_error "unable to restore the default version"

# date & end marker
echo -e "\nFinish: `date`"
echo -e "\n===== End of build and tests ====="
