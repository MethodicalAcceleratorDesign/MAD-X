# run:
# sh scripts/build-test-report.sh [nomail]

# I/O redirection
rm -f build-test-report.log
exec 1> build-test-report.log 2>&1

# env settings
export LC_CTYPE="C"
export PATH=/afs/cern.ch/user/m/mad/madx/madX:$PATH

# setup
thedate=`date "+%Y-%m-%d"`
olddate=`date -d "-50 days" "+%Y-%m-%d"`

# look for failed tests on [lxplus | macosx | windows]
build_test_report ()
{
	if [ ! -s build-test-$1.out ] ; then
		echo "ERROR: report build-test-$1.out not found (or empty) for platform $1"
	else
		cp -f build-test-$1.out tests/reports/${thedate}_build-test-$1.out
		rm -f tests/reports/${olddate}_build-test-$1.out
		perl -ne '/: FAIL|ERROR: / && print' build-test-$1.out > $1-failed.tmp
		if [ -s $1-failed.tmp ] ; then
			perl -ne '/: FAIL|ERROR: / && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-$1.out >> tests-failed.tmp
		fi
		rm -f $1-failed.tmp
	fi
}

# cleaning
rm -f tests-failed.tmp

# lxplus check
build_test_report lxplus

# macosx check
scp -q mad@macserv15865.cern.ch:Projects/madX/build-test-macosx.out .
build_test_report macosx

# report by email if needed
if [ ! -s tests-failed.tmp ] ; then
	echo "ERROR: cannot make build-test-report.out, no input found"
else
	echo "===== Tests Failed ====="                > build-test-report.out
	date                                          >> build-test-report.out
	echo "For details, see report files:"         >> build-test-report.out
	echo "http://cern.ch/madx/madX/tests/reports/${thedate}_build-test-lxplus.out" >> build-test-report.out
	echo "http://cern.ch/madx/madX/tests/reports/${thedate}_build-test-macosx.out" >> build-test-report.out
	echo "http://cern.ch/madx/madX/tests/reports" >> build-test-report.out
	cat tests-failed.tmp                          >> build-test-report.out
	if [ "$1" != "nomail" ] ; then
		cat -v build-test-report.out | mail -s "MAD-X builds and tests report" laurent.deniau@cern.ch # mad-src@cern.ch
		rm -f build-test-report.out
	fi
fi

# cleaning
rm -f tests-failed.tmp

# report errors if any
if [ -s build-test-report.log ] ; then
	cat -v build-test-report.log | mail -s "MAD-X builds and tests errors" laurent.deniau@cern.ch # mad@cern.ch
fi