#!/bin/bash
echo "Executable: $1"
echo "I am at: `pwd`"
echo 'exit;' | $1
