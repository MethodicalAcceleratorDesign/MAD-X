# run:
# sh scripts/build-test-macosx.sh [cleanall]
# tail -f build-test-macosx.out

# I/O redirection
rm -f build-test-macosx.out
exec 1> build-test-macosx.out 2>&1
uname -n > build-test-macosx.run

# env settings
export LC_CTYPE="C"
export PATH=/Users/mad/Projects/madX:/opt/local/bin:$PATH

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
gcc      --version
g++      --version
gfortran --version
make all-macosx-gnu
[ "$?" != "0" ] && echo "ERROR: make all-macosx-gnu failed"

echo -e "\n===== Intel build ====="
icc      --version
ifort    --version
make all-macosx-intel all-macosx
[ "$?" != "0" ] && echo "ERROR: make all-macosx-intel failed"

echo -e "\n===== Binaries dependencies ====="
make infobindep
[ "$?" != "0" ] && echo "ERROR: make infobindep failed"

echo -e "\n===== Tests pointless files ====="
make infotestdep
[ "$?" != "0" ] && echo "ERROR: make infotestdep failed"

echo -e "\n===== Gnu tests (32 bit) ====="
make madx-macosx32-gnu && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-macosx32-gnu failed"

echo -e "\n===== Gnu tests (64 bit) ====="
make madx-macosx64-gnu && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-macosx64-gnu failed"

echo -e "\n===== Intel tests (32 bit) ====="
make madx-macosx32-intel && ls -l madx32 && make cleantest && make tests-all ARCH=32 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-macosx32-intel failed"

echo -e "\n===== Intel tests (64 bit) ====="
make madx-macosx64-intel && ls -l madx64 && make cleantest && make tests-all ARCH=64 NOCOLOR=yes
[ "$?" != "0" ] && echo "ERROR: make tests-all for madx-macosx64-intel failed"

echo -e "\n===== End of build and tests ====="
date

echo "finished" > build-test-macosx.run
