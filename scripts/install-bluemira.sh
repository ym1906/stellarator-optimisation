#!/bin/bash

# read -p "Are you in the correct python environment? (y/n) " answer
# case ${answer:0:1} in
#     y|Y )
#         ;;
#     * )
#         exit;;
# esac

# clean_up() {
#   test -d "$tmp_dir" && rm -rf "$tmp_dir"
# }

# tmp_dir=$( mktemp -d -t install-bm.XXX)
# trap "clean_up $tmp_dir" EXIT

# cd $tmp_dir
INSTALL_CONDA=false
while getopts i: flag
do
    case "${flag}" in
        i) INSTALL_CONDA=true;;
    esac
done

echo
echo Cloning Bluemira...
echo

script_dir=$(dirname "$0")
clone_loc=$script_dir"/../bluemira"
git clone git@github.com:Fusion-Power-Plant-Framework/bluemira.git $clone_loc
cd $clone_loc

echo
echo Getting latest version:
latest_tag=$(git describe --tags $(git rev-list --tags --max-count=1))
echo $latest_tag
echo
git checkout -q $latest_tag

echo
echo Installing...
echo

if [ "$INSTALL_CONDA" = true ] ; then
    source scripts/install-conda.sh
else
    source ~/.miniforge-init.sh
    mamba env create -f conda/environment.yml -n bluemira-stellarator
    mamba activate bluemira-stellarator
fi
pip install -e . --config-settings editable_mode=compat
pre-commit install -f

cd ..
pip install -e .
pre-commit install -f

echo
echo Finished
