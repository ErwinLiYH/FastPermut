# cirrus: module load nvidia/nvhpc-nompi/22.2

# compiler
CC := nvcc
ARCH=sm_70

# compilation flags
CFLAGS := -shared -Xcompiler "-fPIC -Wall -Wextra" -O3 -arch=$(ARCH) -lcublas -lcurand -lcusparse
# the target shared library
TARGET := libpermut_cuda.so
# source files
CUDA_SRCS := permutation_cuda.cu utils_cuda.cu
CPP_SRCS := ../common/utils.cpp
# object files
CUDA_OBJS := $(CUDA_SRCS:.cu=.o)
CPP_OBJS := $(CPP_SRCS:.cpp=.o)
OBJS := $(CUDA_OBJS) $(CPP_OBJS)
# install paths
INSTALL_PATH := $(HOME)
# header files
HEADERS := permutation.cuh

.PHONY: all install uninstall clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cu
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

install:
	mkdir -p $(INSTALL_PATH)/hpcPermut/lib
	cp $(TARGET) $(INSTALL_PATH)/hpcPermut/lib
	mkdir -p $(INSTALL_PATH)/hpcPermut/include
	cp $(HEADERS) $(INSTALL_PATH)/hpcPermut/include

uninstall:
	rm -rf $(INSTALL_PATH)/hpcPermut

clean:
	rm -f $(OBJS) $(TARGET) *.optrpt
