#!/bin/bash

LD_BC_DIR="./OGSTM_BUILD/CMakeFiles/ogstm_lib.dir/src/BC/"
LD_General_DIR="./OGSTM_BUILD/CMakeFiles/ogstm_lib.dir/src/General/"

gfortran -c -std=f2003 -I./OGSTM_BUILD/ -o test_memory_leaks.o test_memory_leaks.f03
gfortran -L${NETCDFF_LIB} -lnetcdff \
    -o test_memory_leaks.x \
    test_memory_leaks.o \
    ${LD_BC_DIR}bc_data.f90.o \
    ${LD_BC_DIR}bc.f90.o \
    ${LD_BC_DIR}rivers.f90.o \
    ${LD_BC_DIR}sponge.f90.o \
    ${LD_BC_DIR}hard_open.f90.o \
    ${LD_BC_DIR}nudging.f90.o \
    ${LD_BC_DIR}bc_aux.o \
    ${LD_General_DIR}calendar.o \
    ${LD_General_DIR}stringop.o

valgrind --leak-check=full ./test_memory_leaks.x

exit 0
