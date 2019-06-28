#!/bin/bash



EXIT_STATUS=0
# 0 = ok
# 1 = not enough arguments
# 2 = wrong argument



# Read parameters

if [ $# -eq 0 ]; then
    echo -e "ERROR: not enough arguments; try $(basename $0) --help"
    EXIT_STATUS=1
    exit $EXIT_STATUS
fi

while [ $# -gt 0 ]; do

    case $1 in

        --non-periodic|-np)
            
            FILE_TIN="files_namelist_riv_default.dat"
            FILE_GIB="files_namelist_gib_default.dat"
            FILE_OPE="files_namelist_dar_default.dat"
            TEST_SUITES="testSuites_default.inc"
            ;;

        --periodic|-p)
            
            FILE_TIN="files_namelist_riv_periodic.dat"
            FILE_GIB="files_namelist_gib_periodic.dat"
            FILE_OPE="files_namelist_dar_periodic.dat"
            TEST_SUITES="testSuites_periodic.inc"
            ;;

        --help|-h)

            echo -e "Usage $(basename $0)"
            echo -e "\t--non-periodic [-np]: links non periodic files"
            echo -e "\t--periodic [-p]: links periodic files"
            echo -e "\t--help [-h]: shows this help"
            exit
            ;;

        *)

            echo -e "ERROR: not a valid option: $1; try $(basename $0) --help"
            EXIT_STATUS=2
            exit $EXIT_STATUS
            ;;

    esac

    shift

done



# Link files accordingly

LINK_TIN="files_namelist_riv.dat"
LINK_GIB="files_namelist_gib.dat"
LINK_OPE="files_namelist_dar.dat"
LINK_SUITES="testSuites.inc"

if [ -L $LINK_TIN ]; then rm $LINK_TIN; fi
if [ -L $LINK_GIB ]; then rm $LINK_GIB; fi
if [ -L $LINK_OPE ]; then rm $LINK_OPE; fi
if [ -L "tests/"$LINK_SUITES ]; then rm "tests/"$LINK_SUITES; fi

ln -sv $FILE_TIN $LINK_TIN
ln -sv $FILE_GIB $LINK_GIB
ln -sv $FILE_OPE $LINK_OPE
ln -sv $TEST_SUITES "tests/"$LINK_SUITES



# Load environment

module load use.own
module load netcdf-ogs pfunit-serial

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

make clean
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
rm "tests/"$LINK_SUITES

module unload netcdf-ogs pfunit-serial
module unload use.own



exit ${EXIT_STATUS}
