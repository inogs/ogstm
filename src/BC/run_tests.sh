#!/bin/bash

EXIT_STATUS=0

make clean # optional

# modify source codes in order to load the right modules
cd src
for f in $(ls *.f03); do
    cat $f | \
        sed 's/bc_aux_mod/bc_aux_testing_mod/g' | \
        sed 's/use modul_param/! use modul_param/g' | \
        sed 's/! integer, parameter :: jp/integer, parameter :: jp/g' \
        > $f".tmp"
    mv $f".tmp" $f
done
cd ..

# actually run tests
make tests
EXIT_STATUS=$?
if [[ ${EXIT_STATUS} -ne 0 ]]; then
    echo -e "ERROR: test script failed with exit status = ${EXIT_STATUS}"
fi

# cleaning up
cd src
for f in $(ls *.f03); do
    cat $f | \
        sed 's/bc_aux_testing_mod/bc_aux_mod/g' | \
        sed 's/! use modul_param/use modul_param/g' | \
        sed 's/integer, parameter :: jp/! integer, parameter :: jp/g' \
        > $f".tmp"
    mv $f".tmp" $f
done
cd ..

exit ${EXIT_STATUS}
