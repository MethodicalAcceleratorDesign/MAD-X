# run:
# sh scripts/build-test-report.sh

export LC_CTYPE="C"

filelist="lxplus-failed.tmp macosx-failed.tmp win-failed.tmp tests-failed.tmp"

# cleaning
rm -f $filelist

# look for failed tests on lxplus
if [ -s build-test-lxplus.out ] ; then
	if [ build-test-lxplus.out -ot build-test-report.out ] ; then
		echo "Warning: build-test-lxplus.out is older than build-test-report.out"
	fi

	perl -ne '/: FAIL/ && print' build-test-lxplus.out > lxplus-failed.tmp

	if [ -s lxplus-failed.tmp ] ; then
		perl -ne '/: FAIL/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-lxplus.out >> tests-failed.tmp
	fi
fi

# look for failed tests on macosx
scp -B mad@macserv15865.cern.ch:/Users/mad/Projects/madX/build-test-macosx.out .

if [ -s build-test-macosx.out ] ; then
	if [ build-test-macosx.out -ot build-test-report.out ] ; then
		echo "Warning: build-test-macosx.out is older than build-test-report.out"
	fi

	perl -ne '/: FAIL/ && print' build-test-macosx.out > macosx-failed.tmp

	if [ -s macosx-failed.tmp ] ; then
		perl -ne '/: FAIL/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-macosx.out >> tests-failed.tmp
	fi
fi

if [ -s tests-failed.tmp ] ; then
	echo "===== Tests Failed =====" >> build-test-report.out
	cat tests-failed.tmp            >> build-test-report.out
	cat -v build-test-report.out | mail -s "MAD-X builds and tests report" mad-src@cern.ch
fi

# cleaning
rm -f $filelist
