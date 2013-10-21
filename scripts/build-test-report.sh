# run:
# sh scripts/build-test-report.sh

# setup
export LC_CTYPE="C"
daynum=`date "+%d"`

# cleaning
rm -f tests-failed.tmp

# look for failed tests on [lxplus | macosx | windows]
build_test_report ()
{
	if [ -s build-test-$1.out ] ; then
		if [ build-test-$1.out -ot build-test-report.out ] ; then
			echo "Warning: build-test-$1.out is older than build-test-report.out"
		fi
		cp -f build-test-$1.out reports/build-test-$1.out.$daynum
		perl -ne '/: FAIL/ && print' build-test-$1.out > $1-failed.tmp
		if [ -s $1-failed.tmp ] ; then
			perl -ne '/: FAIL/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-$1.out >> tests-failed.tmp
		fi
		rm -f $1-failed.tmp
	fi
}

build_test_report lxplus

scp -B mad@macserv15865.cern.ch:/Users/mad/Projects/madX/build-test-macosx.out .
build_test_report macosx

# report by email if needed
if [ -s tests-failed.tmp ] ; then
	echo "===== Tests Failed =====" >> build-test-report.out
	cat tests-failed.tmp            >> build-test-report.out
	cat -v build-test-report.out | mail -s "MAD-X builds and tests report" mad-src@cern.ch
fi

# cleaning
rm -f tests-failed.tmp
