# Include debugging symbols for GDB (Linux, MinGW only) ?
DEBUG_BUILD   := yes

# Compiler command
CC             := g++

# make command
MAKE           := make

# Common CFLAGS
CFLAGS         := -Wall -march=native -mtune=native -pipe

# Use _LINUX for Linux and _WINDOWS for Windows
# Note: Use _WINDOWS for Microsoft Visual Studio
CFLAGS         += -D _LINUX

ifeq ("$(DEBUG_BUILD)","yes")
	CFLAGS += -O0 -g -ggdb
else
	CFLAGS += -O2 -fomit-frame-pointer
endif

# LAPACK
CFLAGS  += $(shell pkg-config --cflags lapack)
LDFLAGS += $(shell pkg-config --libs lapack)

# CFLAGS used for creating compiler generated assembler sources (GCC)
AFLAGS         := $(CFLAGS) -S -masm=intel -fverbose-asm

# Linker flags
LDFLAGS        += -lm -lpthread

# Linker search directories (prepended with -L each)
LDPATH          :=

# Include directories (prepended with -I each)
INC             :=

# Output binary path
BINPATH         := lanczos

# Source directory
SRC             := ../Lanczos

# *.o files directory
OBJ            := obj

# Directory for compiler generated assembler sources
ASM            := asm

CFLAGS += -I"/home/darek/projects/Lanczos/Lanczos"
AFLAGS += -I"/home/darek/projects/Lanczos/Lanczos"

all:
	if [ ! -d "$(OBJ)" ] ; then mkdir -p "$(OBJ)" ; fi
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/file_utils.cpp" -o "$(OBJ)/file_utils.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Lanczno.cpp" -o "$(OBJ)/Lanczno.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/main.cpp" -o "$(OBJ)/main.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Matrixd.cpp" -o "$(OBJ)/Matrixd.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/mrrr.cpp" -o "$(OBJ)/mrrr.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Rescon.cpp" -o "$(OBJ)/Rescon.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/string_utils.cpp" -o "$(OBJ)/string_utils.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/TrMatrixd.cpp" -o "$(OBJ)/TrMatrixd.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/utils.cpp" -o "$(OBJ)/utils.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Vecb.cpp" -o "$(OBJ)/Vecb.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Vecd.cpp" -o "$(OBJ)/Vecd.o"
	$(CC) -c $(CFLAGS) $(INC) "$(SRC)/Veci.cpp" -o "$(OBJ)/Veci.o"
	$(MAKE) link
	$(MAKE) clean

link:
	$(CC) $(CFLAGS) $(INC) \
"$(OBJ)/file_utils.o" \
"$(OBJ)/Lanczno.o" \
"$(OBJ)/main.o" \
"$(OBJ)/Matrixd.o" \
"$(OBJ)/mrrr.o" \
"$(OBJ)/Rescon.o" \
"$(OBJ)/string_utils.o" \
"$(OBJ)/TrMatrixd.o" \
"$(OBJ)/utils.o" \
"$(OBJ)/Vecb.o" \
"$(OBJ)/Vecd.o" \
"$(OBJ)/Veci.o" \
-o "$(BINPATH)" $(LDPATH) $(LDFLAGS)

asm:
	if [ ! -d "$(ASM)" ] ; then mkdir -p "$(ASM)" ; fi
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/file_utils.cpp" -o "$(ASM)/file_utils.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Lanczno.cpp" -o "$(ASM)/Lanczno.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/main.cpp" -o "$(ASM)/main.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Matrixd.cpp" -o "$(ASM)/Matrixd.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/mrrr.cpp" -o "$(ASM)/mrrr.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Rescon.cpp" -o "$(ASM)/Rescon.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/string_utils.cpp" -o "$(ASM)/string_utils.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/TrMatrixd.cpp" -o "$(ASM)/TrMatrixd.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/utils.cpp" -o "$(ASM)/utils.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Vecb.cpp" -o "$(ASM)/Vecb.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Vecd.cpp" -o "$(ASM)/Vecd.asm"
	$(CC) -c $(AFLAGS) $(INC) "$(SRC)/Veci.cpp" -o "$(ASM)/Veci.asm"

clean:
	if [ -w "$(OBJ)/file_utils.o" ] ; then rm "$(OBJ)/file_utils.o" ; fi
	if [ -w "$(OBJ)/Lanczno.o" ] ; then rm "$(OBJ)/Lanczno.o" ; fi
	if [ -w "$(OBJ)/main.o" ] ; then rm "$(OBJ)/main.o" ; fi
	if [ -w "$(OBJ)/Matrixd.o" ] ; then rm "$(OBJ)/Matrixd.o" ; fi
	if [ -w "$(OBJ)/mrrr.o" ] ; then rm "$(OBJ)/mrrr.o" ; fi
	if [ -w "$(OBJ)/Rescon.o" ] ; then rm "$(OBJ)/Rescon.o" ; fi
	if [ -w "$(OBJ)/string_utils.o" ] ; then rm "$(OBJ)/string_utils.o" ; fi
	if [ -w "$(OBJ)/TrMatrixd.o" ] ; then rm "$(OBJ)/TrMatrixd.o" ; fi
	if [ -w "$(OBJ)/utils.o" ] ; then rm "$(OBJ)/utils.o" ; fi
	if [ -w "$(OBJ)/Vecb.o" ] ; then rm "$(OBJ)/Vecb.o" ; fi
	if [ -w "$(OBJ)/Vecd.o" ] ; then rm "$(OBJ)/Vecd.o" ; fi
	if [ -w "$(OBJ)/Veci.o" ] ; then rm "$(OBJ)/Veci.o" ; fi

.PHONY: asm
