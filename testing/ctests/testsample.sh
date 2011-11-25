#!/bin/bash

echo "Executable: $1"

echo "ctests directory: $2"

./$1 < $2/job.sample.madx
