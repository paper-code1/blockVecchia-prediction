#----------------complier configurations-------------------------
_SUPPORT_VECCHIA_?=TRUE
# _USE_MAGMA_?=TRUE
#specify cuda directory
# _CUDA_ROOT_=$(CUDA_HOME)
# _CUDA_ARCH_ ?= 70

#----------------compliers -------------------------
CXX = g++ # 10.2.0
# NVCC=$(_CUDA_ROOT_)/bin/nvcc

# NVOPTS = -ccbin $(CXX) --compiler-options -fno-strict-aliasing
CXXFLAGS = -O2 -Wall -std=c++17 -fopenmp -Wsign-compare

ifdef _DEBUG_
  CXXFLAGS += -g -Xcompiler -rdynamic
  # NVOPTS += -G -g -lineinfo
else
  CXXFLAGS += -O3
  # NVOPTS += -O3
endif

# ifdef _USE_MAGMA_
  # CXXFLAGS += -DUSE_MAGMA
  # _MAGMA_ROOT_?=$(HOME)/dev/magma-2.7.2
  # NVOPTS += -DUSE_MAGMA
# endif

#---------------- directory -------------------------
# object and bin files
OBJ_DIR=./obj
BIN_DIR=./bin
SRC_DIR=./src
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


INCLUDES=
INCLUDES+= -I.
INCLUDES+= -I./include
# INCLUDES+= -I${CUDA_ROOT}/include
# INCLUDES+= -I${NLOPT_ROOT}/include
INCLUDES+= -I${GSL_ROOT}/include # used for matern kernel, bessel function

LIB_PATH=
# LIB_PATH+= -L${CUDA_ROOT}/lib64
# LIB_PATH+= -L${NLOPT_ROOT}/lib
LIB_PATH+= -L${GSL_ROOT}/lib

# ifdef _USE_MAGMA_
# 	LIB_PATH+= -L${_MAGMA_ROOT_}/lib
# 	INCLUDES+= -I$(_MAGMA_ROOT_)/include
# endif

# libraries to link against
LIB= -lm 
# LIB+= -lnlopt  
LIB+= -lgsl #2.6
LIB+= -lgslcblas
LIB+= -lblas
# ifdef _USE_MAGMA_
# 	LIB+= -lmagma -lcusparse
# endif
# LIB+= -lcublas -lcudart
LIB+= -lgomp # 4.1.0
LIB+= -lstdc++ # 10.2.0

#---------------- make -------------------------


TARGET = $(BIN_DIR)/predict
# Source files
SRC = main.cpp kmeans.cpp csv_utils.cpp cluster_utils.cpp
# Object files
OBJS = $(addprefix $(OBJ_DIR)/, $(SRC:.cpp=.o))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIB) $(LIB_PATH)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET)

.PHONY: clean
