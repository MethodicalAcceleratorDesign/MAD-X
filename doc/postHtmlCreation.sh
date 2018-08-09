CURR=$(pwd)
OPHYS=$CURR/html/userguide
cd $OPHYS
cp manual_temp.html manual.html
for FILE in *temp.html; do
	echo " cleaning html $FILE"
	#sed -i 's/\.<\/span><\/td\>/\\<\/span\>\<\/td\>/g' $FILE
	sed -i 's/\.<\/span><\/td>/\<\/span><\/td>/g' $FILE
	
done
