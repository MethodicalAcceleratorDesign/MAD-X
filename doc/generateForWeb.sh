#!/bin/bash

# Generates a web based on the latex file. 
# Requires the latexml package to be installed.
# Written by Veronica Berglyd Olsen, Feb 2018 and adapted to the need for MAD-X by Tobias Persson Aug 2018
#

CURR=$(pwd)

TUSER=$CURR/user_manual_temp
MPHYS=$CURR/latexuguide


OPHYS=$CURR/html/webguide
OBUILD=$CURR/html/build_full

# LaTeXML Options
MATHJAX='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML'
FORMAT=html5


#mkdir -pv $OUSERS
mkdir -pv $OPHYS

cd $MPHYS/figures
./generatePlots.sh 
cd ../..

echo ""
echo "*******************************"
echo "* Generating User Manual HTML *"
echo "*******************************"
echo ""

# Make temp directory
rm -rfv $TUSER
rsync -avPh $MPHYS/ $TUSER
cd $TUSER


# Remove unsupported stuff for LaTeX source
for FILE in *.tex; do
    echo "Cleaning up file $FILE"
    sed -i 's/\\begin{cverbatim}/\\begin{verbatim}/g' $FILE
    sed -i 's/\\end{cverbatim}/\\end{verbatim}/g' $FILE
    sed -i 's/\\begin{ctverbatim}/\\begin{verbatim}/g' $FILE
    sed -i 's/\\end{ctverbatim}/\\end{verbatim}/g' $FILE
    sed -i 's/\\begin{longtabu}/\\begin{tabular}/g' $FILE
    sed -i 's/\\end{longtabu}/\\end{tabular}/g' $FILE
    sed -i 's/\\include{coverpage}/\\input{coverpage.tex}/g' $FILE
    sed -i 's/\\\[/$/g' $FILE
    sed -i 's/\\\]/$/g' $FILE   
    sed -i 's/\\\end{minipage}}/\\\end{minipage}/g' $FILE
    sed -i 's/\madbox/madboxweb/g' $FILE
  #  sed -i 's/\/]/$/g' $FILE

#    sed -i 's/\\ttitem/\\item/g' $FILE
    sed -i 's/\\arraybackslash//g' $FILE
    sed -i '/\\todo/d' $FILE
    sed -i '/\\pdfbookmark/d' $FILE
    sed -i '/\\pdfbook/d' $FILE
    sed -i '/\\listoftables/d' $FILE
    sed -i '/\\listoffigure/d' $FILE
    sed -i '/\\framebox\[\\textwidth\]/d' $FILE
    sed -i '/\\usepackage\[totoc\]{idxlayout}/d' $FILE
    sed -i '/\\usepackage{makeidx}/d' $FILE
    sed -i '/\\printindex/d' $FILE


done



# Build


echo ""
echo "**********************************"
echo "* Generating Physics Manual HTML *"
echo "**********************************"
echo ""

cd $TUSER
#make
#cp $MPHYS/sixphys.pdf $CURR/html/physics_manual.pdf
latexml uguide.tex --includestyles | latexmlpost --dest=$OPHYS/manual_temp.html --format=$FORMAT --javascript=$MATHJAX -
#$CURR/cleanupHTML.py $OPHYS
#rm -v $OPHYS/*.html
#echo "<?php header('Location: manual.php'); ?>" > $OPHYS/index.php
./postHtmlCreation.sh

echo ""
echo "**********"
echo "*  DONE  *"
echo "**********"
echo ""
echo "The content of the folder 'html' can now be uploaded to "
echo ""
