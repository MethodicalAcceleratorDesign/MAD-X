# run:
# sh scripts/build-test-report.sh

list="lxplus-failed.tmp macosx-failed.tmp win-failed.tmp tests-failed.tmp build-test-failed.tmp"

rm -f $list

# look for failed tests on lxplus
if [ -s build-test-lxplus.out ] ; then
	perl -ne '/: \W\[31mFAIL\W\[0m/ && print' build-test-lxplus.out > lxplus-failed.tmp
	if [ -s lxplus-failed.tmp ] ; then
		perl -ne '/: \W\[31mFAIL\W\[0m/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-lxplus.out >> tests-failed.tmp
	fi
fi

# look for failed tests on macosx
if [ -s build-test-macosx.out ] ; then
	rm -f macosx-failed.tmp
	perl -ne '/: \W\[31mFAIL\W\[0m/ && print' build-test-macosx.out > macosx-failed.tmp
	if [ -s macosx-failed.tmp ] ; then
		perl -ne '/: \W\[31mFAIL\W\[0m/ && print ; /-> (madx-\S+)/ && print "\n$1:\n"' build-test-macosx.out >> tests-failed.tmp
	fi
fi

if [ -s tests-failed.tmp ] ; then
	echo "===== Tests Failed =====" >> build-test-failed.tmp
	cat tests-failed.tmp            >> build-test-failed.tmp
	mail -s "MAD-X build and tests report" laurent.deniau@cern.ch < build-test-failed.tmp
fi

# cat build-test-failed.tmp

# cleaning
rm -f $list

