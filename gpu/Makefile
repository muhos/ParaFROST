################################################################################
#
# Muhammad Osama Mahmoud, Eindhoven University of Technology (TU/e).
# Mathematics and Computer Science Department, SET cluster.
# GEARS Project, All rights reserved.
#
################################################################################

# CUDA path
CUDA_PATH ?= /usr/local/cuda-10.1

# operating system
HOST_OS   := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
TARGET_OS ?= $(HOST_OS)
ifeq (,$(filter $(TARGET_OS),linux))
    $(error ERROR - unsupported value $(TARGET_OS) for TARGET_OS!)
endif

# target arch
TARGET_ARCH   := $(shell uname -m)

# host compiler
HOST_COMPILER := g++

# device compiler
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

# target size
TARGET_SIZE   := $(shell getconf LONG_BIT)

# compiler flags
NVCCFLAGS := -m${TARGET_SIZE}
CCFLAGS   := -std=c++14 -Wno-format-zero-length
MTHREADS  := -pthreads 
LDFLAGS   := 

# Debug build flags
WARN := -w 
ifeq ($(dbg),1)
      NVCCFLAGS += -g
      BUILD_TYPE := debug
      WARN :=
else  ifeq ($(assert),1)
      NVCCFLAGS += -O3 -use_fast_math
      BUILD_TYPE := release
      WARN :=
else
      NVCCFLAGS += -O3 -use_fast_math -DNDEBUG
      BUILD_TYPE := release
endif

# combine all flags
ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))

# Common includes 
INCLUDES  := -I../../common/inc
LIBRARIES :=

# Gencode arguments
SMS ?= 35 50 60 70 75
All_GOOD := 1
ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified <<<)
All_GOOD := 0
endif

ifeq ($(GENCODE_FLAGS),)
  # Generate SASS code for each SM architecture in $(SMS)
  $(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
  # Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
  HIGHEST_SM := $(lastword $(sort $(SMS)))
  ifneq ($(HIGHEST_SM),)
  GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
  endif
endif

ifeq ($(All_GOOD),0)
EXEC ?= @echo "[@]"
endif

# Target rules
all: build

build: parafrost

check.deps:
ifeq ($(All_GOOD),0)
	@echo "Compilation will be skipped due to the above missing dependencies"
else
	@echo "Compilation is ready - all dependencies have been met"
endif

pfalloc.o:pfalloc.cu
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<
 
pfkernels.o:pfkernels.cu
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<
 
pfsimp.o:pfsimp.cu
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

pfsolve.o:pfsolve.cpp
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

pfpdm.o:pfpdm.cpp
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

pfargs.o:pfargs.cpp
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

pfmain.o:pfmain.cpp
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

parafrost: pfargs.o pfalloc.o pfkernels.o pfsimp.o pfpdm.o pfsolve.o pfmain.o
	$(EXEC) $(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -o $@ $+ $(LIBRARIES)
	$(EXEC) mkdir -p build/$(BUILD_TYPE)
	$(EXEC) cp $@ build/$(BUILD_TYPE)

clean:
	rm -f pfalloc.o pfkernels.o pfargs.o pfsimp.o pfpdm.o pfsolve.o pfmain.o parafrost  
	rm -rf build/$(BUILD_TYPE)/parafrost
