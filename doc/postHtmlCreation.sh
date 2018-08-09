#!/bin/bash
CURR=$(pwd)
OPHYS=$CURR/html/userguide
cd $OPHYS
cp manual.html manual_temp.html
for FILE in *temp.html; do
	echo " cleaning html $FILE"
	#sed -i 's/\.<\/span><\/td\>/\\<\/span\>\<\/td\>/g' $FILE
	sed -i 's/\.<\/span><\/td>/\<\/span><\/td>/g' $FILE
	echo "\.\<\/span\>\<\/td\>"
done