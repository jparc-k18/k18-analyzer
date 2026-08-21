##
# @file		Makefile
# @author	S.H. Hayakawa
# @date		2015/03/10
# @brief	Makefile for K1.8 offline analyzer
#
##### configuration ##################################################
CXX	= g++
#
BIN_DIR = bin
SRC_DIR	= src
INC_DIR	= include
USR_DIR	= usr
LIB_DIR	= lib
SRC_BLD_DIR	= $(SRC_DIR)/build
USR_BLD_DIR	= $(USR_DIR)/build
#
# ROOT
root_config	= root-config
ROOT_INCLUDE	= $(shell $(root_config) --cflags)
ROOT_LIBS	= $(shell $(root_config) --libs)
#
# HD DAQ Unpacker
unpacker_install	?= $(HOME)/work/k18analyzer/e05/unpacker
unpacker_include_dir	:= $(unpacker_install)/include
unpacker_include	:= -I$(unpacker_include_dir)
unpacker_lib_dir	:= $(unpacker_install)/lib
unpacker_libs		:= -Wl,-rpath,$(unpacker_lib_dir) \
	-L$(unpacker_lib_dir) -lHDDAQUnpacker -lUnpackerManager
#
CFLAGS		= -g -O
# -Wall
#DFLAGS		= -Df2cFortran -Dextname -DDEBUG
FLAGS	= $(CFLAGS) $(DFLAGS) -I. $(ROOT_INCLUDE) \
		$(unpacker_include) $(addprefix -I,$(SRC_DIR) $(INC_DIR))

LIBOBJ	= $(LIB_DIR)/libSksAna.a
SKSLIB	= -L$(LIB_DIR)/ -lSksAna
LIBS	= $(ROOT_LIBS) $(unpacker_libs) -lUnpackerManager
LIBSONL	= $(ROOT_LIBS) $(unpacker_libs)

SKSSRCS	= $(wildcard $(SRC_DIR)/*.cc)
HEADS	= $(wildcard $(SRC_DIR)/*.hh)
USRSRCS	= $(wildcard $(USR_DIR)/*.cc)
DST_DIR	= dst
DST_BLD_DIR	= $(DST_DIR)/build
DSTSRCS	= $(wildcard $(DST_DIR)/*.cc)
DSTDEPS	= $(DSTSRCS:$(DST_DIR)/%.cc=$(DST_BLD_DIR)/%.d)
DSTOBJ	= $(DSTSRCS:$(DST_DIR)/%.cc=$(DST_BLD_DIR)/%.o)
DSTTARGET	= $(DSTSRCS:$(DST_DIR)/%.cc=$(BIN_DIR)/%)

SKSDEPS	= $(SKSSRCS:$(SRC_DIR)/%.cc=$(SRC_BLD_DIR)/%.d)
USRDEPS	= $(USRSRCS:$(USR_DIR)/%.cc=$(USR_BLD_DIR)/%.d)

SKSOBJ	= $(SKSSRCS:$(SRC_DIR)/%.cc=$(SRC_BLD_DIR)/%.o)
USROBJ	= $(USRSRCS:$(USR_DIR)/%.cc=$(USR_BLD_DIR)/%.o)
TARGET	= $(USRSRCS:$(USR_DIR)/User%.cc=$(BIN_DIR)/%)

##### target #########################################################
.PHONY: all lib usr clean

all: lib usr dst

-include $(SKSDEPS) $(USRDEPS) $(DSTDEPS)

usr: $(TARGET)
dst: $(DSTTARGET)
lib: $(LIBOBJ)

$(LIBOBJ): $(SKSOBJ)
	@echo -e "\e[33;1m=== Archiving $@ \e[m"
	@mkdir -p $(LIB_DIR)
	ar r $(LIBOBJ) $(SKSOBJ)

$(TARGET): $(LIBOBJ) $(USROBJ)
	@echo -e "\e[35;1m=== Linking $@ \e[m"
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $(FLAGS) $(@:$(BIN_DIR)/%=$(USR_BLD_DIR)/User%.o) \
	$(SKSLIB) $(LIBS)

$(SRC_BLD_DIR)/%.o: $(SRC_DIR)/%.cc
	@echo -e "\e[32;1m=== Compiling $@ \e[m"
	@mkdir -p $(SRC_BLD_DIR)
	$(CXX) $(FLAGS) -o $@ -MMD -c $<

$(USR_BLD_DIR)/%.o: $(USR_DIR)/%.cc
	@echo -e "\e[32;1m=== Compiling $@ \e[m"
	@mkdir -p $(USR_BLD_DIR)
	$(CXX) $(FLAGS) -o $@ -MMD -c $<

$(DSTTARGET): $(LIBOBJ) $(DSTOBJ)
	@echo -e "\e[35;1m=== Linking $@ \e[m"
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $(FLAGS) $(@:$(BIN_DIR)/%=$(DST_BLD_DIR)/%.o) \
	$(SKSLIB) $(LIBS)

$(DST_BLD_DIR)/%.o: $(DST_DIR)/%.cc
	@echo -e "\e[32;1m=== Compiling $@ \e[m"
	@mkdir -p $(DST_BLD_DIR)
	$(CXX) $(FLAGS) -o $@ -MMD -c $<

clean:
	@echo -e "\e[32;1m=== Cleaning \e[m"
	@rm -rfv $(USROBJ) $(LIBOBJ)
	@find . \( -name "*~" -o -name "\#*\#" \) -exec rm -fv \{\} \;

distclean:
	@echo -e "\e[32;1m=== Cleaning \e[m"
	@rm -rfv $(TARGET) $(SRC_BLD_DIR) $(USR_BLD_DIR)
	@find . \( -name "*~" -o -name "\#*\#" \) -exec rm -fv \{\} \;
