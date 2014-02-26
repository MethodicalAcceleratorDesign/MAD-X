# run:
# sh scripts/test-valgrind.sh [madx]

readonly testdir="tests-valgrind"
readonly summary="valgrind-summary.txt"
readonly pattern="== Invalid (read|write) of size"

# select madx, check for existence
madx="$1"
[ "$madx" == "" ] && madx="madx" 
[ ! -x "$madx" ] && echo "error: invalid program '$madx' for tests with valgrind" && exit 1

# go to tests directory
cd $testdir
[ "$?" != "0" ] && echo "error: copy 'tests' tree into '$testdir' before running tests with valgrind" && exit 1
echo "moved to `pwd`"

# clean summary, set the date
echo "tests started at `date`" > $summary

# run all tests with valgrind
for i in test-*; do
  echo "running test $i (produce $i.valgrind): `date`"
  cd $i
  valgrind -v --leak-check=full --track-origins=yes ../../$madx $i.madx > $i.valgrind 2>&1
  grep -E "$pattern" $i.valgrind /dev/null >> ../$summary
  cd ..
done

cd ..