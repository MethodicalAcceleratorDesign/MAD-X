#!/bin/bash
echo "Binary postfix: $1"
echo 'exit;' | ./madx$1
