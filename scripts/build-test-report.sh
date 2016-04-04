#! /bin/bash
# run:
# bash scripts/build-test-report.sh [noecho] [clean|cleanall] [forcereport] [force] [nomail]

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

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
    if [ "$1" == "clean" -o "$1" == "cleanall" ] ; then
        clean="clean"
    fi

# force to build the report even if remote processes are not finished or errors occur
    if [ "$1" == "force"  ] ; then
        force="force"
    fi

# force to (re)build the report even if non empty build-test-report.out exists
    if [ "$1" == "forcereport"  ] ; then
        forcereport="forcereport"
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

# clean all and die
[ "$clean" == "clean" ] && rm -f build-test-*.tmp build-test-*.out && die

# prevent erasing existing report even with "force", override with "forcereport" (for debugging)
[ "$forcereport" != "forcereport" -a -s "build-test-report.out" ] && exit

# general settings
readonly thedate=`date "+%Y-%m-%d"`
readonly webdir="http://cern.ch/madx/madX"

readonly windir="mad@macserv15865W10.cern.ch:madX"
readonly linuxdir="mad@macserv15865LX.cern.ch:madX"
readonly macosxdir="mad@macserv15865.dyndns.cern.ch:madX"

# error handler
check_error ()
{
    [ "$?" != "0" ] && echo "ERROR: $1" && die
}

# clean tempory files
clean_tmp ()
{
    rm -f build-test-*.tmp
}

# clear reports older than 51 days
clear_old_reports ()
{
    find tests/reports -ctime +51  -name '*_build-test-[lmw][aix]*.out' -exec rm {} \;
    find tests/reports -ctime +501 -name '*_build-test-report.out' -exec rm {} \;
}

# check for completed jobs [lxplus | macosx | linux | win]
build_test_completed ()
{
    local marker="not found"

    for arch in "$@" ; do
        if [ -s build-test-$arch.out ] ; then
            marker=`perl -ne '/===== End of build and tests =====/ && print "found"' build-test-$arch.out`
            check_error "unable to search for end marker (perl)"
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
    local dir
    for arch in "$@" ; do
        eval dir=\$${arch}dir
        scp -q -p "$dir/build-test-$arch.out" build-test-$arch.out
        check_error "unable to retrieve $arch remote report (scp)"
        if [ -s build-test-$arch.out ] ; then
            cat build-test-$arch.out | tr -d '\r' > build-test-$arch.tr
            mv -f build-test-$arch.tr build-test-$arch.out
            # retrieve binaries for download of last builds
            scp -q -p "$dir/madx-${arch}64-gnu*" .
            scp -q -p "$dir/madx-${arch}32-gnu*" .
            scp -q -p "$dir/numdiff-${arch}64-gnu*" .
            scp -q -p "$dir/numdiff-${arch}32-gnu*" .
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
        echo -e "\n=====\n$webdir/tests/reports/${thedate}_build-test-$arch.out$completed" >> build-test-result.tmp

        if [ -s build-test-$arch.out ] ; then
            cp -f build-test-$arch.out tests/reports/${thedate}_build-test-$arch.out
            check_error "backup of build-test-$arch.out failed (cp)"

            perl -ne '/: FAIL$|^ERROR: / && print' build-test-$arch.out >> build-test-failed.tmp
            check_error "unable to search for failures or errors (perl)"

            perl -ne '/: FAIL$|^ERROR: / && print ;
                      /===== Testing (madx-\S+) =====/ && print "\n$1:\n"' build-test-$arch.out >> build-test-result.tmp
            check_error "unable to build report summary (perl)"
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

    echo "===== Tests $status ====="                                           > build-test-report.out
    date                                                                      >> build-test-report.out
    svn info svn+ssh://svn.cern.ch/reps/madx/trunk/madX                       >> build-test-report.out
    echo "For details, see report files:"                                     >> build-test-report.out
    echo "$webdir/tests"                                                      >> build-test-report.out
    echo "$webdir/tests/reports"                                              >> build-test-report.out
    echo "$webdir/tests/reports/${thedate}_build-test-report.out (this file)" >> build-test-report.out
    cat build-test-result.tmp                                                 >> build-test-report.out

    cp -f build-test-report.out tests/reports/${thedate}_build-test-report.out
    check_error "backup of build-test-report.out failed (cp)"

    if [ "$nomail" != "nomail" ] ; then
        cat -v build-test-report.out | mail -s "MAD-X builds and tests report ${thedate}: $status" mad-src@cern.ch
        check_error "unable to email report summary (mail)"
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
    check_error "unable to email report errors (check mail)"
fi

clear_old_reports
clean_tmp
