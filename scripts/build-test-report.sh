#! /bin/bash
# run:
# bash scripts/build-test-report.sh [noecho] [clean] [force] [doreport] [nomail]

# env settings
export PATH="`pwd`:$PATH"

# I/O redirection
rm -f build-test-report.log
if [ "$1" = "noecho" ] ; then
  shift
  exec &> build-test-report.log
else
  exec 2>&1 | tee build-test-report.log
fi

# parse arguments
while [ "$1" != "" ] ; do

# clean tmp files and remote reports, and exit if not 'force'd to continue
  if [ "$1" == "clean" ] ; then
    clean="clean"
  fi

# force to build the report even if remote processes are not finished or errors occur
  if [ "$1" == "force"  ] ; then
    force="force"
  fi

# force to (re)build the report even if non empty build-test-report.out exists (for debugging)
  if [ "$1" == "doreport"  ] ; then
    doreport="doreport"
  fi

# do not send the report by email (for debugging)
  if [ "$1" == "nomail"  ] ; then
    nomail="nomail"
  fi

  shift
done

# quit on error if not forced to continue
die ()
{
  [ "$force" != "force" ] && exit 1
}

# error handler
check_error ()
{
  if [ "$?" != "0" ] ; then
    echo -e "\nERROR: $1"
    [ "$2" != "no-exit" ] && die
  fi
}

# go to madx reports directory
[ -d madx-reports ] && cd madx-reports && echo "moving down to madx-reports"

# clean all and die
[ "$clean" == "clean" ] && rm -f build-test-*.tmp build-test-*.out && die

# prevent erasing existing report even with "force", override with "doreport"
[ "$doreport" != "doreport" -a -s "build-test-report.out" ] && exit

# general settings
readonly repdir="mad@lxplus.cern.ch:madx/madx-reports"
readonly thedate=`date "+%Y-%m-%d"`

readonly winsrc="mad@macserv15865W10.cern.ch:"
readonly linuxsrc="mad@macserv15865LX.cern.ch:"
readonly macosxsrc="mad@macserv15865.cern.ch:"
readonly lxplussrc="mad@lxplus.cern.ch:madx"

# clean tempory files
clean_tmp ()
{
  rm -f build-test-*.tmp
}

# keep 1 year of reports
clear_old_reports ()
{
  find db -ctime +365 -name '*_build-test-[lmw][aix]*.out' -exec rm {} \;
  find db -ctime +365 -name '*_build-test-report.out' -exec rm {} \;
}

# check for completed jobs [lxplus | macosx | linux | win]
build_test_completed ()
{
  local marker="not found"

  for arch in "$@" ; do
    if [ -s build-test-$arch.out ] ; then
      marker=`perl -ne '/===== End of build and tests =====/ && print "found"' build-test-$arch.out`
      check_error "unable to search for end marker (perl)" "no-exit"
    fi
    [ "$marker" != "found" ] && incomplete="incomplete, " && return 1
  done
  return 0
}

build_test_check ()
{
  [ "$force" != "force" ] && build_test_completed "$@" || die
  return 0
}

# retrieve remote report [lxplus | macosx | linux | win]
build_test_remote ()
{
  local src
  for arch in "$@" ; do
    eval src=\$${arch}src
    scp -q -p "$src/build-test-$arch.out" build-test-$arch.out
    check_error "unable to retrieve $arch remote report (scp)" "no-exit"
    if [ -s build-test-$arch.out ] ; then
      cat build-test-$arch.out | tr -d '\r' > build-test-$arch.tr
      mv -f build-test-$arch.tr build-test-$arch.out
      # remove local copies to ensure proper scp (no -force option)
      rm -f madx-${arch}64-gnu* madx-${arch}32-gnu*
      rm -f numdiff-${arch}64-gnu* numdiff-${arch}32-gnu*
      # retrieve binaries for download of last builds
      scp -q -p "$src/madx-${arch}64-gnu*" .
      scp -q -p "$src/madx-${arch}32-gnu*" .
      scp -q -p "$src/numdiff-${arch}64-gnu*" .
      scp -q -p "$src/numdiff-${arch}32-gnu*" .
    fi
  done
  return 0
}

# look for failed tests [lxplus | macosx | linux | win]
build_test_report ()
{
  local completed

  for arch in "$@" ; do
    build_test_completed $arch && completed="" || completed=" (not found or incomplete)"
    echo -e "\n=====\n$repdir/db/${thedate}_build-test-$arch.out$completed" >> build-test-result.tmp

    if [ -s build-test-$arch.out ] ; then
      cp -f build-test-$arch.out db/${thedate}_build-test-$arch.out
      check_error "backup of build-test-$arch.out failed (cp)" "no-exit"

      perl -ne '/: FAIL$|^ERROR: / && print' build-test-$arch.out >> build-test-failed.tmp
      check_error "unable to search for failures or errors (perl)" "no-exit"

      perl -ne '/: FAIL$|^ERROR: / && print ;
           /===== Testing (madx-\S+) =====/ && print "\n$1:\n"' build-test-$arch.out >> build-test-result.tmp
      check_error "unable to build report summary (perl)" "no-exit"
    else
      echo "ERROR: report build-test-$arch.out not found (or empty) for platform $arch"
    fi
  done
  return 0
}

# send daily reports summary by email [lxplus | macosx | linux | win]
build_test_send ()
{
  local status
  local count

  if [ -s build-test-failed.tmp ] ; then
    count=`wc -l < build-test-failed.tmp`
    status="${incomplete}failed ($count)"
  elif [ "$incomplete" != "" ] ; then
    status="${incomplete}failed"
  else
    status="passed"
  fi

  git fetch

  echo "===== Tests $status ====="                                > build-test-report.out
  date                                                           >> build-test-report.out
  git show --summary origin/master                               >> build-test-report.out
#    svn info svn+ssh://svn.cern.ch/reps/madx/trunk/madX         >> build-test-report.out
  echo "For details, see report files:"                          >> build-test-report.out
  echo "$repdir/db/${thedate}_build-test-report.out (this file)" >> build-test-report.out
  cat build-test-result.tmp                                      >> build-test-report.out

  cp -f build-test-report.out db/${thedate}_build-test-report.out
  check_error "backup of build-test-report.out failed (cp)" "no-exit"

  if [ "$nomail" != "nomail" ] ; then
    cat -v build-test-report.out | mail -s "MAD-X builds and tests report ${thedate}: $status" mad-src@cern.ch
    check_error "unable to email report summary (mail)" "no-exit"
  fi
  return 0
}

# cleaning
clean_tmp

# check if local reports are finished
build_test_check  lxplus

# retrieve remote reports
build_test_remote        macosx linux win

# check if non-local reports are finished
build_test_check         macosx linux win

# build the final report
build_test_report lxplus macosx linux win

# send the final report
build_test_send   lxplus macosx linux win

# report errors by email if any
if [ "$nomail" != "nomail" -a -s build-test-report.log ] ; then
  cat -v build-test-report.log | mail -s "MAD-X builds and tests report errors (${thedate})" mad@cern.ch
  check_error "unable to email report errors (check mail)" "no-exit"
fi

# backup last-build
if [ -x madx-linux64-gnu -a -x ../releases ] ; then
 cp -f madx-linux64-gnu ../releases
fi

clear_old_reports
clean_tmp
