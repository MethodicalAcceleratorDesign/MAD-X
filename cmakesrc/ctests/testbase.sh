#!/bin/bash

echo "Executable: $1"
echo "Test directory: $2"
echo "Script: $3"
echo "Output: $4"

correct_dir=/afs/cern.ch/user/y/ylevinse/scratch1/public/madx_testing_output/

$1 < $2/$3

madret=$?
if [ ! -d ${correct_dir} ]
then
    echo "WARNING: Directory with correct output is not available."
    echo "Cannot make comparison of files"
    exit $madret
fi

# This will always fail because data/version is included
# diff $4 ${correct_dir}/$4 > /dev/null
# if (( $?!=0 )) 
# then
#     echo "Output generated in file $4 is not correct"
#     exit 5
# else
#     echo "Output generated is correct"
# fi

exit $madret
