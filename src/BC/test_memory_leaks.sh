#!/bin/bash

INCLUDE_DIR="./src/"
LD_BC_DIR="./src/"
LD_General_DIR="./src/"

gfortran -c -std=f2003 -I${INCLUDE_DIR} -o test_memory_leaks.o test_memory_leaks.f03
gfortran -L${NETCDFF_LIB} -lnetcdff \
    -o test_memory_leaks.x \
    test_memory_leaks.o \
    ${LD_BC_DIR}bc_data.o \
    ${LD_BC_DIR}bc.o \
    ${LD_BC_DIR}rivers.o \
    ${LD_BC_DIR}sponge.o \
    ${LD_BC_DIR}hard_open.o \
    ${LD_BC_DIR}nudging.o \
    ${LD_BC_DIR}bc_aux_testing.o \
    ${LD_General_DIR}calendar.o \
    ${LD_General_DIR}stringop.o

valgrind --leak-check=full --show-leak-kinds=all ./test_memory_leaks.x

exit 0
