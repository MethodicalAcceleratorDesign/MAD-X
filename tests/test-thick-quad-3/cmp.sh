for PARAM1 in 5 50 500 5000 ; do
echo "NSLICES = $PARAM1"
sed "s/PARAM1/$PARAM1/" fivecell.madx > tmp
time ../../madx < tmp > out_$PARAM1
done

echo "THICK"
time ../../madx test-thick-quad-3.madx > out_thick
