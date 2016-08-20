#!/bin/bash
set -eu

if [ $# -gt 0 ] ; then
	OUF="$1"
else
	OUF="Makefile"
fi


OUF=`realpath "${OUF}"`

cat Makefile.in > "${OUF}"

echo "Output file: \"${OUF}\""

SRC_DIR=`grep "SRC[[:blank:]]*:=" Makefile.in | awk -F '=' '{print $2}' | tr -d ' '`
SRC_DIR=`realpath "${SRC_DIR}"`
echo "Source directory: \"${SRC_DIR}\""

echo "CFLAGS += -I\"${SRC_DIR}\"" >> "${OUF}"
echo "AFLAGS += -I\"${SRC_DIR}\"" >> "${OUF}"
echo "" >> "${OUF}"

cd "${SRC_DIR}"

CPP_FILES=`ls *.cpp`

echo -e "all:" >> "${OUF}"
echo -e "\tif [ ! -d \"\$(OBJ)\" ] ; then mkdir -p \"\$(OBJ)\" ; fi" >> "${OUF}"

for f in $CPP_FILES ; do
	OBJ_FILE=`echo "${f}" | sed 's/.cpp/.o/'`
	echo -e "\t\$(CC) -c \$(CFLAGS) \$(INC) \"\$(SRC)/${f}\" -o \"\$(OBJ)/${OBJ_FILE}\"" >> "${OUF}"
done

echo -e "\t\$(MAKE) link" >> "${OUF}"
echo -e "\t\$(MAKE) clean" >> "${OUF}"
echo    "" >> "${OUF}"
echo    "link:" >> "${OUF}"
echo -e "\t\$(CC) \$(CFLAGS) \$(INC) \\" >> "${OUF}"

for f in $CPP_FILES ; do
	OBJ_FILE=`echo "${f}" | sed 's/.cpp/.o/'`
	echo    "\"\$(OBJ)/${OBJ_FILE}\" \\" >> "${OUF}"
done

echo    "-o \"\$(BINPATH)\" \$(LDPATH) \$(LDFLAGS)" >> "${OUF}"
echo    "" >> "${OUF}"
echo    "asm:" >> "${OUF}"
echo -e "\tif [ ! -d \"\$(ASM)\" ] ; then mkdir -p \"\$(ASM)\" ; fi" >> "${OUF}"

for f in $CPP_FILES ; do
	ASM_FILE=`echo "${f}" | sed 's/.cpp/.asm/'`
	echo -e "\t\$(CC) -c \$(AFLAGS) \$(INC) \"\$(SRC)/${f}\" -o \"\$(ASM)/${ASM_FILE}\"" >> "${OUF}"
done

echo    "" >> "${OUF}"
echo    "clean:" >> "${OUF}" >> "${OUF}"

for f in $CPP_FILES ; do
	OBJ_FILE=`echo "${f}" | sed 's/.cpp/.o/'`
	echo -e "\tif [ -w \"\$(OBJ)/${OBJ_FILE}\" ] ; then rm \"\$(OBJ)/${OBJ_FILE}\" ; fi" >> "${OUF}"
done

echo    "" >> "${OUF}"
echo    ".PHONY: asm" >> "${OUF}"

exit 0
