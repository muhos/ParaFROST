################################################################################
#
# Muhammad Osama Mahmoud, Eindhoven University of Technology (TU/e).
# Mathematics and Computer Science Department, MDSE, SET Group.
# GEARS Project, All rights reserved.
#
################################################################################

# operating system
HOST_OS   := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
TARGET_OS ?= $(HOST_OS)
ifeq (,$(filter $(TARGET_OS),linux))
    $(error ERROR - unsupported value $(TARGET_OS) for TARGET_OS!)
endif

# target arch
TARGET_ARCH := $(shell uname -m)

# host compiler
HOST_COMPILER ?= g++

# target size
TARGET_SIZE := $(shell getconf LONG_BIT)

# compiler  flags
MTHREADS  := -pthreads
CCFLAGS   := -m${TARGET_SIZE} -std=c++14 
LDFLAGS   := 

# Debug build flags
WARN := -w 
ifeq ($(dbg),1)
      CCFLAGS += -g
      BUILD_TYPE := debug
      WARN :=
else  ifeq ($(assert),1)
      CCFLAGS += -O3 -use_fast_math
      BUILD_TYPE := release
      WARN :=
else
      CCFLAGS += -O3 -use_fast_math -DNDEBUG
      BUILD_TYPE := release
endif

# Common includes 
INCLUDES  := 
LIBRARIES :=

# Target rules
all: build

build: parafrost

sigma.o:sigma.cpp
	$(EXEC) $(HOST_COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $(WARN) $<

solver.o:solver.cpp
	$(EXEC) $(HOST_COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $(WARN) $<

pdm.o:pdm.cpp
	$(EXEC) $(HOST_COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $(WARN) $<

args.o:args.cpp
	$(EXEC) $(HOST_COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $(WARN) $<

main.o:main.cpp
	$(EXEC) $(HOST_COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $(WARN) $<

parafrost: args.o pdm.o sigma.o solver.o main.o
	$(EXEC) $(HOST_COMPILER) $(LDFLAGS) $(CCFLAGS) -o $@ $+ $(LIBRARIES)
	$(EXEC) mkdir -p $(TARGET_ARCH)/$(TARGET_OS)/$(BUILD_TYPE)
	$(EXEC) cp $@ $(TARGET_ARCH)/$(TARGET_OS)/$(BUILD_TYPE)

clean:
	rm -f args.o sigma.o solver.o pdm.o main.o parafrost  
	rm -rf $(TARGET_ARCH)/$(TARGET_OS)/$(BUILD_TYPE)/parafrost

clobber: clean
