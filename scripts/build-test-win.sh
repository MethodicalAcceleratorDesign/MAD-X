#! /bin/bash
# run:
# bash scripts/build-test-win.sh [noecho] [clone|update] [cleanall] [notest]

# env settings
export PATH="`pwd`:/c/Program Files/gnuplot/bin:$PATH"

# error handler
check_error ()
{
	if [ "$?" != "0" ] ; then
		echo -e "\nERROR: $1"
    [ "$2" != "no-exit" ] && exit 1
	fi
}

# change directory
if [ "$PWD" != "$HOME/madx-nightly" ] ; then
  cd "$HOME/madx-nightly"
  check_error "unable to move to madx-nightly directory"
  export PATH="`pwd`:$PATH"
fi

# set env 32 or 64 bit
set_env ()
{
    local lo="s/mingw(32|64)/mingw$1/g"
    local up="s/MINGW(32|64)/MINGW$1/g"

    export PATH=`echo "$PATH"                       | sed -r -e "$lo"`
    export MANPATH=`echo $MANPATH                   | sed -r -e "$lo"`
    export MSYSTEM=`echo "$MSYSTEM"                 | sed -r -e "$up"`
    export ACLOCAL_PATH=`echo "$ACLOCAL_PATH"       | sed -r -e "$lo"`
    export PKG_CONFIG_PATH=`echo "$PKG_CONFIG_PATH" | sed -r -e "$lo"`
}

# I/O redirection
rm -f build-test-win.out
if [ "$1" = "noecho" ] ; then
	shift
  export NOCOLOR=yes
	exec > build-test-win.out 2>&1
	check_error "redirection with noecho failed"
else
	exec > >(tee build-test-win.out) 2>&1
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

echo -e "\n===== Git clone/update ====="
if [ "$1" = "clone" ] ; then
  shift
  update="clone"
  rm -rf madx-nightly && \
  git clone https://github.com/MethodicalAcceleratorDesign/MAD-X.git madx-nightly && \
  cd madx-nightly
  check_error "git clone failed"
elif [ "$1" = "update" ] ; then
  shift # faster "clone"
  update="update"
  [ -d madx-nightly ] && cd madx-nightly && echo "moving down to madx-nightly"
  git fetch --tags && \
  git reset --hard origin/master && \
  git clean -xdf
  check_error "git update failed"
else
  update=""
  echo "Skipped (no explicit request)."
fi

echo -e "\n===== Git info ====="
git log -1 --format="Branch:  %d%nCommit:   %H%nAuthor:   %an <%ae>%nDate:     %ad (%ar)%nSubject:  %s"

lasttest=`git for-each-ref refs/tags --sort=-committerdate --format='%(refname)' --count=1`
echo -e "\nFiles changed since (last) release: ${lasttest##*/}"
git diff --name-status $lasttest

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
set_env 64
gcc      --version
g++      --version
gfortran --version

set_env 32
make all-win32-gnu
check_error "make all-win32-gnu failed" "no-exit"

set_env 64
make all-win64-gnu
check_error "make all-win64-gnu failed" "no-exit"

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

	echo -e "\n===== Testing madx-win32-gnu ====="
	make madx-win32-gnu && ls -l madx32 && make cleantest && make tests-all COMP=gnu ARCH=32 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-win32-gnu failed" "no-exit"

	echo -e "\n===== Testing madx-win64-gnu ====="
	make madx-win64-gnu && ls -l madx64 && make cleantest && make tests-all COMP=gnu ARCH=64 NOCOLOR=$NOCOLOR
	check_error "make tests-all for madx-win64-gnu failed" "no-exit"
fi

# restore the default version
make madx-win32-gnu > /dev/null && make madx-win64-gnu > /dev/null
check_error "unable to restore the default version"

# date & end marker
echo -e "\nFinish: `date`"
echo -e "\n===== End of build and tests ====="
