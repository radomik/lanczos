#!/bin/sh

root="../"
src_files=`find "${root}" \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.py" -o -name "Makefile" -o -name "*.sh" -o -name "*.makefile" -o -name "Makefile*" -o -name "*.txt" -o -name "*.m" -o -name "*.c" \) | tr -s '\n' ' '`
geany -i -s ${src_files}

