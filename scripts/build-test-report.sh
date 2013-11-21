# run:
# sh scripts/build-test-report.sh

# I/O redirection
rm -f build-test-report.log
exec 1> build-test-report.log 2>&1

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

# setup
readonly thedate=`date "+%Y-%m-%d"`
readonly olddate=`date -d "-50 days" "+%Y-%m-%d"`
readonly nomail="$1"
readonly macdir="mad@macserv15865.cern.ch:Projects/madX"
readonly lxpdir="http://cern.ch/madx/madX"

clean_tmp ()
{
	rm -f *.tmp
}

clean_exit ()
{
	clean_tmp
	exit
}

# check for finished jobs [lxplus | macosx | win]
build_test_check ()
{
	for arch in "$@" ; do
		if [ ! -s build-test-$arch.run ] || [ "`cat build-test-$arch.run`" != "finished" ] ; then
			clean_exit
		fi
	done
}

# look for failed tests [lxplus | macosx | win]
build_test_report ()
{
	for arch in "$@" ; do
		if [ ! -s build-test-$arch.out ] ; then
			echo "ERROR: report build-test-$arch.out not found (or empty) for platform $arch"
		else
			rm -f tests/reports/${olddate}_build-test-$arch.out
			cp -f build-test-$arch.out tests/reports/${thedate}_build-test-$arch.out
			[ "$?" != "0" ] && echo "ERROR: backup of build-test-$arch.out failed (check cp)"

			perl -ne '/: FAIL|ERROR: / && print' build-test-$arch.out >> tests-failed.tmp
			[ "$?" != "0" ] && echo "ERROR: unable to search for failures or errors (check perl)"

			perl -ne '/: FAIL|ERROR: / && print ; /===== Testing (madx-\S+) =====/ && print "\n$1:\n"' build-test-$arch.out >> tests-result.tmp
			[ "$?" != "0" ] && echo "ERROR: unable to build report summary (check perl)"
		fi
	done
}

# send daily reports summary by email [lxplus | macosx | win]
build_test_send ()
{
	local status

	if [ -s tests-failed.tmp ] ; then
		status="failed"
	else
		status="passed"
	fi

	echo "===== Tests $status ====="                                  > build-test-report.out
	date                                                             >> build-test-report.out
	echo "For details, see report files:"                            >> build-test-report.out
	echo "$lxpdir/tests/reports/${thedate}_build-test-report.out"    >> build-test-report.out
	for arch in "$@" ; do
		echo "$lxpdir/tests/reports/${thedate}_build-test-$arch.out" >> build-test-report.out
	done
	echo "$lxpdir/tests/reports"                                     >> build-test-report.out
	echo "$lxpdir/tests"                                             >> build-test-report.out
	cat tests-result.tmp                                             >> build-test-report.out

	if [ "$nomail" != "nomail" ] ; then
		cat -v build-test-report.out | mail -s "MAD-X builds and tests report ($status)" mad-src@cern.ch
		[ "$?" != "0" ] && echo "ERROR: unable to email report summary (check mail)"
	fi
	cp -f build-test-report.out tests/reports/${thedate}_build-test-report.out
}

# tag reports as processed [lxplus | macosx | win]
build_test_proc ()
{
	for arch in "$@" ; do
		echo "processed `date`" > build-test-$arch.run
	done
}

# cleaning
clean_tmp

# check if local reports are finished
build_test_check  lxplus

# retrieve remote reports
scp -q "$macdir/build-test-macosx.*" "$macdir/build-test-win.*" .
[ "$?" != "0" ] && echo "ERROR: unable to retrieve macosx report (check scp)"

# check if non-local reports are finished
build_test_check         macosx win

# build the final report
build_test_report lxplus macosx win

# send the final report
build_test_send   lxplus macosx win

# mark all reports as processed
build_test_proc   lxplus macosx win

# update status of remote reports
scp -q build-test-macosx.run build-test-win.run "$macdir"
[ "$?" != "0" ] && echo "ERROR: unable to update macosx report (check scp)"

# report errors if any
if [ -s build-test-report.log ] ; then
	cat -v build-test-report.log | mail -s "MAD-X builds and tests report errors" mad@cern.ch
	[ "$?" != "0" ] && echo "ERROR: unable to email report errors (check mail)"
fi

clean_exit
