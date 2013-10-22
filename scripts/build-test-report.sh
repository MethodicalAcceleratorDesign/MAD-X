# run:
# sh scripts/build-test-report.sh

# setup
export LC_CTYPE="C"
thedate=`date "+%Y-%m-%d"`

# look for failed tests on [lxplus | macosx | windows]
build_test_report ()
{
	if [ -s build-test-$1.out ] ; then
		cp -f build-test-$1.out tests/reports/${thedate}_build-test.out.bak
		perl -ne '/: FAIL/ && print' build-test-$1.out > $1-failed.tmp
		if [ -s $1-failed.tmp ] ; then
			perl -ne '/: FAIL/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-$1.out >> tests-failed.tmp
		fi
		rm -f $1-failed.tmp
	fi
}

# cleaning
rm -f tests-failed.tmp

# lxplus check
build_test_report lxplus

# macosx check
scp -q -B mad@macserv15865.cern.ch:/Users/mad/Projects/madX/build-test-macosx.out .
build_test_report macosx

# report by email if needed
if [ -s tests-failed.tmp ] ; then
	echo "===== Tests Failed =====" >  build-test-report.out
	date                            >> build-test-report.out
	cat tests-failed.tmp            >> build-test-report.out
	if [ "$1" != "nomail" ] ; then
		cat -v build-test-report.out | mail -s "MAD-X builds and tests report" mad-src@cern.ch
		rm -f build-test-report.out
	fi
fi

# cleaning
rm -f tests-failed.tmp
