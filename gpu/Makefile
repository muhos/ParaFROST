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
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

# target size
TARGET_SIZE   := $(shell getconf LONG_BIT)

# compiler flags
CCSTD     := -std=c++14
NVCCFLAGS := -m${TARGET_SIZE}
CCFLAGS   := $(CCSTD) -Wno-format-zero-length
LDFLAGS   := 

ifneq ($(MAKECMDGOALS),clean)
      $(info found operating system "$(HOST_OS)")
      $(info building with "nvcc (master) + $(HOST_COMPILER) (slave)" compiler)
      $(info for a "$(TARGET_ARCH) $(TARGET_SIZE)-bit" architecture)
endif

# Debug build flags
WARN := -w 
ifeq ($(debug),1)
      NVCCFLAGS += -g
      BUILD_TYPE := debug
      WARN :=
      ifneq ($(MAKECMDGOALS),clean)
        $(info using "$(NVCCFLAGS) $(CCFLAGS)" for a debug version with assertions)
      endif
else  ifeq ($(assert),1)
      NVCCFLAGS += -O3 -use_fast_math
      BUILD_TYPE := release
      WARN :=
      ifneq ($(MAKECMDGOALS),clean)
        $(info using "$(NVCCFLAGS) $(CCFLAGS)" for a release version with assertions)
      endif
else
      NVCCFLAGS += -O3 -use_fast_math -DNDEBUG
      BUILD_TYPE := release
      ifneq ($(MAKECMDGOALS),clean)
        $(info using "$(NVCCFLAGS) $(CCFLAGS)" for a release version)
      endif
endif

# combine all flags
ALL_CCFLAGS := $(WARN)
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
BIN       := parafrost
LIB       := lib$(BIN).a
CPPOBJEXT := o
CUOBJEXT  := cuo 


mainsrc   := pfmain
cusrc     := $(sort $(wildcard *.cu))
allcppsrc := $(sort $(wildcard *.cpp))
cppsrc    := $(filter-out $(mainsrc).cpp,$(allcppsrc))
cuobj     := $(patsubst %.cu,%.$(CUOBJEXT),$(cusrc))
cppobj    := $(patsubst %.cpp,%.$(CPPOBJEXT),$(cppsrc))

ifneq ($(MAKECMDGOALS),clean)
	ifeq ($(cusrc),)
		$(error No CUDA source files exist)
	endif
	ifeq ($(cppsrc),)
		$(error No C++ source files exist)
	endif
	ifeq ($(cuobj),)
		$(error No CUDA object files to generate)
	endif
	ifeq ($(cppobj),)
		$(error No C++ object files to generate)
	endif
    $(info compiling all sources to generate "$(BIN)" and its library "$(LIB)")
endif

all: $(BIN)

$(LIB): $(cuobj) $(cppobj)
	ar rc $@ $+
	ranlib $@

$(BIN): $(mainsrc).$(CPPOBJEXT) $(LIB)
	$(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -o $@ $(mainsrc).$(CPPOBJEXT) -L. -l$(BIN) $(LIBRARIES)
	mkdir -p build/$(BUILD_TYPE)
	cp $@ build/$(BUILD_TYPE)

$(mainsrc).$(CPPOBJEXT): $(mainsrc).cpp
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $< 

%.$(CUOBJEXT): %.cu
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $< 
	
%.$(CPPOBJEXT): %.cpp
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

clean:
	rm -f *.$(CPPOBJEXT) *.$(CUOBJEXT) $(LIB) $(BIN)  
	rm -rf build/$(BUILD_TYPE)/$(BIN)
