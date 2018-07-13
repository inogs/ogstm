#!/bin/bash

gfortran -c -std=f2003 -I./src/ -o test_memory_leaks.o test_memory_leaks.f03
gfortran -L./src/ -o test_memory_leaks.x test_memory_leaks.o src/bc_data.o src/calendar.o src/stringop.o
valgrind --leak-check=full ./test_memory_leaks.x

exit 0
