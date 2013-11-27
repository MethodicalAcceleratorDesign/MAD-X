# run:
# sh scripts/build-test-report.sh

# I/O redirection
rm -f build-test-report.log
exec 1> build-test-report.log 2>&1

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

# exclusive settings
readonly nomail="$1"
readonly force="$1"

# general settings
readonly thedate=`date "+%Y-%m-%d"`
readonly olddate=`date -d "-50 days" "+%Y-%m-%d"`
readonly macdir="mad@macserv15865.cern.ch:Projects/madX"
readonly lxpdir="http://cern.ch/madx/madX"

clean_tmp ()
{
	rm -f build-test-*.tmp
}

clean_all ()
{
	rm -f build-test-*.tmp build-test-*.out
}

# check for completed jobs [lxplus | macosx | win]
build_test_completed ()
{
	local marker

	for arch in "$@" ; do
		if [ -s build-test-$arch.out ] ; then
			marker=`perl -ne '/===== End of build and tests =====/ && print "found"' build-test-$arch.out`
			[ "$?" != "0" ] && echo "ERROR: unable to search for end marker (perl)"
			[ "$marker" != "found" ] && return 1
		fi
	done
	return 0
}

build_test_check ()
{
	if [ "$force" != "force" ] ; then
		if ! build_test_completed "$@" ; then
			echo "Reports not yet completed, retrying later..."
			echo `date`
			clean_tmp && exit
		fi
	fi
}

# look for failed tests [lxplus | macosx | win]
build_test_report ()
{
	local completed

	for arch in "$@" ; do
		build_test_completed $arch && completed="" || completed=" (incomplete)"
		echo "$lxpdir/tests/reports/${thedate}_build-test-$arch.out$completed" >> build-test-result.tmp


		if [ ! -s build-test-$arch.out ] ; then
			echo "ERROR: report build-test-$arch.out not found (or empty) for platform $arch"
		else
			rm -f tests/reports/${olddate}_build-test-$arch.out
			cp -f build-test-$arch.out tests/reports/${thedate}_build-test-$arch.out
			[ "$?" != "0" ] && echo "ERROR: backup of build-test-$arch.out failed (cp)"

			perl -ne '/: FAIL|ERROR: / && print' build-test-$arch.out >> build-test-failed.tmp
			[ "$?" != "0" ] && echo "ERROR: unable to search for failures or errors (perl)"

			perl -ne '/: FAIL|ERROR: /                 && print ;
			          /===== Testing (madx-\S+) =====/ && print "\n$1:\n"' build-test-$arch.out >> build-test-result.tmp
			[ "$?" != "0" ] && echo "ERROR: unable to build report summary (perl)"
		fi
	done
}

# send daily reports summary by email [lxplus | macosx | win]
build_test_send ()
{
	local status

	[ -s build-test-failed.tmp ] && status="failed" || status="passed"

	echo "===== Tests $status ====="                               > build-test-report.out
	date                                                          >> build-test-report.out
	echo "For details, see report files:"                         >> build-test-report.out
	echo "$lxpdir/tests/reports/${thedate}_build-test-report.out" >> build-test-report.out
	echo "$lxpdir/tests/reports"                                  >> build-test-report.out
	echo "$lxpdir/tests"                                          >> build-test-report.out
	cat build-test-result.tmp                                     >> build-test-report.out

	if [ "$nomail" != "nomail" ] ; then
		cat -v build-test-report.out | mail -s "MAD-X builds and tests report ${thedate}: $status" mad-src@cern.ch
		[ "$?" != "0" ] && echo "ERROR: unable to email report summary (mail)"
	fi
	cp -f build-test-report.out tests/reports/${thedate}_build-test-report.out
}

# clean all and quit
[ "$1" == "clean" -o "$1" == "cleanall" ] && clean_all && exit

# report already exists, prevent erasing it even with "force"
[ -s "build-test-report.out" ] && exit

# cleaning
clean_tmp

# check if local reports are finished
build_test_check  lxplus

# retrieve remote reports
scp -q -p "$macdir/build-test-macosx.out" "$macdir/build-test-win.out" .
[ "$?" != "0" ] && echo "ERROR: unable to retrieve macosx report (scp)"

# check if non-local reports are finished
build_test_check         macosx win

# build the final report
build_test_report lxplus macosx win

# send the final report
build_test_send   lxplus macosx win

# report errors if any
if [ "$nomail" != "nomail" -a -s build-test-report.log ] ; then
	cat -v build-test-report.log | mail -s "MAD-X builds and tests report errors (${thedate})" mad@cern.ch
	[ "$?" != "0" ] && echo "ERROR: unable to email report errors (check mail)"
fi

clean_tmp
