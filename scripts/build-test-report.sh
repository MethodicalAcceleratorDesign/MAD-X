# run:
# sh scripts/build-test-report.sh [clean] [forcereport] [force] [nomail]

# I/O redirection
rm -f build-test-report.log
exec 1> build-test-report.log 2>&1

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

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
readonly olddate=`date -d "-50 days" "+%Y-%m-%d"` 	# linux
# readonly olddate=`date -v-50d "+%Y-%m-%d"` 		# macosx
readonly srvdir="mad@macserv15865.cern.ch:Projects/madX"
readonly webdir="http://cern.ch/madx/madX"

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

# check for completed jobs [lxplus | macosx | win]
build_test_completed ()
{
	local marker="not found"

	for arch in "$@" ; do
		if [ -s build-test-$arch.out ] ; then
			marker=`perl -ne '/===== End of build and tests =====/ && print "found"' build-test-$arch.out`
			check_error "unable to search for end marker (perl)"
		fi
		[ "$marker" != "found" ] && return 1
	done
	return 0
}

build_test_check ()
{
	[ "$force" != "force" ] && build_test_completed "$@" || die
	return 0
}

# retrieve remote report [lxplus | macosx | win]
build_test_remote ()
{
	for arch in "$@" ; do
		scp -q -p "$srvdir/build-test-$arch.out" build-test-$arch.out
		check_error "unable to retrieve $arch remote report (scp)"
		if [ ! -s build-test-$arch.out ] ; then
			cat build-test-$arch.out | tr -d '\r' > build-test-$arch.out.tr
			mv -f build-test-$arch.out.tr build-test-$arch.out
		fi
	done
	return 0
}

# look for failed tests [lxplus | macosx | win]
build_test_report ()
{
	local completed

	for arch in "$@" ; do
		build_test_completed $arch && completed="" || completed=" (not found or incomplete)"
		echo -e "\n=====\n$webdir/tests/reports/${thedate}_build-test-$arch.out$completed" >> build-test-result.tmp

		if [ ! -s build-test-$arch.out ] ; then
			echo "ERROR: report build-test-$arch.out not found (or empty) for platform $arch"
		else
			rm -f tests/reports/${olddate}_build-test-$arch.out
			cp -f build-test-$arch.out tests/reports/${thedate}_build-test-$arch.out
			check_error "backup of build-test-$arch.out failed (cp)"

			perl -ne '/: FAIL|ERROR:|error: / && print' build-test-$arch.out >> build-test-failed.tmp
			check_error "unable to search for failures or errors (perl)"

			perl -ne '/: FAIL|ERROR:|error: / && print ;
			          /===== Testing (madx-\S+) =====/ && print "\n$1:\n"' build-test-$arch.out >> build-test-result.tmp
			check_error "unable to build report summary (perl)"
		fi
	done
	return 0
}

# send daily reports summary by email [lxplus | macosx | win]
build_test_send ()
{
	local status
	local count

	if [ -s build-test-failed.tmp ] ; then
		count=`wc -l < build-test-failed.tmp`
		status="failed ($count)"
	else
		status="passed"
	fi

	echo "===== Tests $status ====="                                           > build-test-report.out
	date                                                                      >> build-test-report.out
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
build_test_remote 		 macosx win

# check if non-local reports are finished
build_test_check         macosx win

# build the final report
build_test_report lxplus macosx win

# send the final report
build_test_send   lxplus macosx win

# report errors by email if any
if [ "$nomail" != "nomail" -a -s build-test-report.log ] ; then
	cat -v build-test-report.log | mail -s "MAD-X builds and tests report errors (${thedate})" mad@cern.ch
	check_error "unable to email report errors (check mail)"
fi

clean_tmp
