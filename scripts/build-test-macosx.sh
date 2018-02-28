#! /bin/bash
# run:
# bash scripts/build-test-macosx.sh [noecho] [clone|update|clean] [nobuild] [notest]

# env settings
export PATH="`pwd`:/opt/local/bin:$PATH"

# error handler
check_error ()
{
  if [ "$?" != "0" ] ; then
    echo -e "\nERROR: $1"
    [ "$2" != "no-exit" ] && exit 1
  fi
}

# I/O redirection
rm -f build-test-macosx.out
if [ "$1" = "noecho" ] ; then
  shift
  export NOCOLOR=yes
  exec > build-test-macosx.out 2>&1
  check_error "redirection with noecho failed"
else
  exec > >(tee build-test-macosx.out) 2>&1
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
echo "PWD   : `echo $PWD`"

echo -e "\n===== Git clone/update/clean ====="
if [ "$1" = "clone" ] ; then
  shift # git clone
  rm -rf madx-nightly && \
  git clone https://github.com/MethodicalAcceleratorDesign/MAD-X.git madx-nightly
  check_error "git clone failed"
  [ -d madx-nightly ] && cd madx-nightly && echo "moving down to cloned madx-nightly"

elif [ "$1" = "update" ] ; then
  shift # faster "clone" + git cleanup
  [ -d madx-nightly ] && cd madx-nightly && echo "moving down to updated madx-nightly"
  git fetch --tags && \
  git reset --hard origin/master
  check_error "git update failed"
  git clean -dfqx
  check_error "git cleanup failed" "no-exit"

elif [ "$1" = "clean" ] ; then
  shift # git cleanup
  [ -d madx-nightly ] && cd madx-nightly && echo "moving down to cleaned madx-nightly"
  git clean -fqx
  check_error "git cleanup failed" "no-exit"

else
  echo "Skipped (no explicit request)."
fi

echo -e "\n===== Git info ====="
git log -1 --format="Branch:  %d%nCommit:   %H%nAuthor:   %an <%ae>%nDate:     %ad (%ar)%nSubject:  %s"

lasttest=`git for-each-ref refs/tags --sort=-committerdate --format='%(refname)' --count=1`
echo -e "\nFiles changed since (last) release: ${lasttest##*/}"
git diff --name-status $lasttest

################################################################################
if [ "$1" = "nobuild" ] ; then
  echo -e "\nBuild and tests skipped (explicit request)."
  echo -e "\nFinish: `date`"
  echo -e "\n===== End ====="
  exit
fi
################################################################################

echo -e "\n===== Release number ====="
cat VERSION

echo -e "\n===== Gnu build ====="
gcc      --version
g++      --version
gfortran --version

make all-macosx32-gnu
check_error "make all-macosx32-gnu failed" "no-exit"

make all-macosx64-gnu
check_error "make all-macosx64-gnu failed" "no-exit"

echo -e "\n===== Intel build ====="

if [ "`which icc`" != "" -a "`which ifort`" != "" ] ; then
  icc      --version
  ifort    --version

  make all-macosx32-intel
  check_error "make all-macosx32-intel failed" "no-exit"

  make all-macosx64-intel
  check_error "make all-macosx64-intel failed" "no-exit"
else
  echo "Intel compilers not found, skipped."
fi

echo -e "\n===== Binaries dependencies ====="
make infobindep
check_error "make infobindep failed" "no-exit"

echo -e "\n===== Tests pointless files ====="
make cleantest && make infotestdep
check_error "make infotestdep failed" "no-exit"

echo -e "\n===== Running tests (long) ====="
if [ "$1" = "notest" ] ; then
  shift
  echo "Skipped (explicit request)."
else
  echo ""

if [ "`which icc`" != "" -a "`which ifort`" != "" ] ; then
  echo -e "\n===== Testing madx-macosx32-intel ====="
  make madx-macosx32-intel && ls -l madx32 && make cleantest && make tests-all COMP=intel ARCH=32 NOCOLOR=$NOCOLOR
  check_error "make tests-all for madx-macosx32-intel failed" "no-exit"

  echo -e "\n===== Testing madx-macosx64-intel ====="
  make madx-macosx64-intel && ls -l madx64 && make cleantest && make tests-all COMP=intel ARCH=64 NOCOLOR=$NOCOLOR
  check_error "make tests-all for madx-macosx64-intel failed" "no-exit"
else
  echo "Intel compilers not found, skipped."
fi

  echo -e "\n===== Testing madx-macosx32-gnu ====="
  make madx-macosx32-gnu && ls -l madx32 && make cleantest && make tests-all COMP=gnu ARCH=32 NOCOLOR=$NOCOLOR
  check_error "make tests-all for madx-macosx32-gnu failed" "no-exit"

  echo -e "\n===== Testing madx-macosx64-gnu ====="
  make madx-macosx64-gnu && ls -l madx64 && make cleantest && make tests-all COMP=gnu ARCH=64 NOCOLOR=$NOCOLOR
  check_error "make tests-all for madx-macosx64-gnu failed" "no-exit"
fi

# restore the default version
make madx-macosx32-gnu > /dev/null && make madx-macosx64-gnu > /dev/null
check_error "unable to restore the default version" "no-exit"

# date & end marker
echo -e "\nFinish: `date`"
echo -e "\n===== End of build and tests ====="

exit
