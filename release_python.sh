#!/bin/bash
set -euo pipefail; IFS=$'\n\t'

NAME=$( python setup.py --name )
VER=$( python setup.py --version )

echo "========================================================================"
echo "Tagging $NAME v$VER"
echo "========================================================================"

#git tag v$VER
#git push origin v$VER

echo "========================================================================"
echo "Releasing $NAME v$VER on PyPI"
echo "========================================================================"

python setup.py sdist
twine upload dist/*
rm -r dist/ *.egg-info
