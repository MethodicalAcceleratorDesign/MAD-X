#!/bin/bash

echo "Binary postfix: $1"

echo "ctests directory: $2"

./madx$1 < $2/job.sample.madx
