CURR=$(pwd)
OPHYS=$CURR/html/webguide
cd $OPHYS
cp manual_temp.html manual.html
for FILE in manual.html; do
	echo " cleaning html $FILE"
	#sed -i 's/\.<\/span><\/td\>/\\<\/span\>\<\/td\>/g' $FILE
	sed -i 's/\.<\/span><\/td>/\<\/span><\/td>/g' $FILE
	sed -i 's/Untitled Document/MAD-X guide/g' $FILE
    sed -i 's/Untitled Document/MAD-X guide/g' $FILE
	sed -i  s'/http:\/\/madx.web.cern.ch\/madx\/madX\/doc\/usrguide\///g' $FILE
done
