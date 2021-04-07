################################################################################
#
# Muhammad Osama Mahmoud, Eindhoven University of Technology (TU/e).
# Mathematics and Computer Science Department, SET cluster.
# GEARS Project, All rights reserved.
#
################################################################################

# CUDA path
CUDA_PATH ?= $(CUDA_HOME)

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
NVCC          := $(CUDA_HOME)/bin/nvcc -ccbin $(HOST_COMPILER)

# target size
TARGET_SIZE   := $(shell getconf LONG_BIT)

# compiler flags
NVCCFLAGS := -m${TARGET_SIZE}
CCFLAGS   := -std=c++14 -Wno-format-zero-length
MTHREADS  := -pthreads 
LDFLAGS   := 

# Debug build flags
WARN := -w 
ifeq ($(debug),1)
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
SMS ?= 60 70 75
ifeq ($(SMS),)
$(error No SM architectures have been specified)
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

# Target rules
cusrc  := $(sort $(wildcard *.cu))
cppsrc := $(sort $(wildcard *.cpp))
obj    := $(patsubst %.cu,%.o,$(cusrc))
obj    += $(patsubst %.cpp,%.o,$(cppsrc))

ifeq ($(cusrc),)
$(error No CUDA source files exist)
endif

ifeq ($(cppsrc),)
$(error No C++ source files exist)
endif

ifeq ($(obj),)
$(error No object files to generate)
endif

all: build

build: parafrost

libpfrost.a: $(obj)
	ar rc $@ $+
	ranlib $@

parafrost: pfmain.o libpfrost.a
	$(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -o $@ pfmain.o -L. -lpfrost $(LIBRARIES)
	mkdir -p build/$(BUILD_TYPE)
	cp $@ build/$(BUILD_TYPE)

pfmain.o:pfmain.cpp
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

%.o:%.cu
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<
	
%.o:%.cpp
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $(WARN) $<

clean:
	rm -f *.o libpfrost.a parafrost  
	rm -rf build/$(BUILD_TYPE)/parafrost
