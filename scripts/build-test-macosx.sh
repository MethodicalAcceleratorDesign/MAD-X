# run:
# sh scripts/build-test-macosx.sh [cleanall]
# tail -f build-test-macosx.out

# env settings
export LC_CTYPE="C"
export PATH="/Users/mad/Projects/madX:/opt/local/bin:$PATH"

# I/O redirection
rm -f build-test-macosx.out
exec 1> build-test-macosx.out 2>&1

# error handler
check_error ()
{
	if [ "$?" != "0" ] ; then
		echo "ERROR: $1"
		exit 1
	fi
}

echo "\n===== Start of build and tests ====="
date
uname -m -n -r -s

echo "\n===== SVN update ====="
svn update
check_error "svn update failed"

echo "\n===== Release number ====="
cat VERSION

echo "\n===== Clean build ====="
if [ "$1" = "cleanall" ] ; then
	make cleanall
	check_error "make cleanall failed"
else
	echo "Skipped (no explicit request)."
fi 

echo "\n===== Gnu build ====="
gcc      --version
g++      --version
gfortran --version
make all-macosx-gnu
# to handle bad fortran compiler, restart from scratch (only once)
if [ "$?" != "0" ] ; then
	make cleanall && make cleanall ARCH=32
	check_error "make cleanall failed"
	make all-macosx-gnu
	check_error "make all-macosx-gnu failed"
fi

echo "\n===== Intel build ====="
icc      --version
ifort    --version
make all-macosx-intel all-macosx
check_error "make all-macosx-intel failed"

echo "\n===== Binaries dependencies ====="
make infobindep
check_error "make infobindep failed"

echo "\n===== Tests pointless files ====="
make cleantest && make infotestdep
check_error "make infotestdep failed"

echo "\n===== Testing madx-macosx64-intel ====="
make madx-macosx64-intel && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
check_error "make tests-all for madx-macosx64-intel failed"

echo "\n===== Testing madx-macosx32-intel ====="
make madx-macosx32-intel && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
check_error "make tests-all for madx-macosx32-intel failed"

echo "\n===== Testing madx-macosx64-gnu ====="
make madx-macosx64-gnu && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
check_error "make tests-all for madx-macosx64-gnu failed"

echo "\n===== Testing madx-macosx32-gnu ====="
make madx-macosx32-gnu && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
check_error "make tests-all for madx-macosx32-gnu failed"

# restore the default version
make madx-macosx32 > /dev/null && make madx-macosx64 > /dev/null
check_error "unable to restore the default version"

# date & end marker
date
echo "\n===== End of build and tests ====="
