# cirrus: module load nvidia/nvhpc-nompi/22.2

# compiler
CC := nvcc
ARCH=sm_70

# compilation flags
CFLAGS := -shared -Xcompiler "-fPIC -Wall -Wextra" -O3 -arch=$(ARCH) -lcublas -lcurand -lcusparse -L../common -lpermut_utils
# the target shared library
TARGET := libpermut_cuda.so
# source files
SRCS := permutation_cuda.cu utils_cuda.cu
# object files
OBJS := $(SRCS:.cu=.o)
# install paths
ifndef PERMUT_PATH
	INSTALL_PATH := $(HOME)
else
	INSTALL_PATH := $(PERMUT_PATH)
endif
# header files
HEADERS := permutation.cuh

.PHONY: all install uninstall clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cu
	$(CC) $(CFLAGS) -c -o $@ $<

install:
	mkdir -p $(INSTALL_PATH)/hpcPermut/lib
	cp $(TARGET) $(INSTALL_PATH)/hpcPermut/lib
	mkdir -p $(INSTALL_PATH)/hpcPermut/include/cuda
	cp $(HEADERS) $(INSTALL_PATH)/hpcPermut/include/cuda

uninstall:
	rm -rf $(INSTALL_PATH)/hpcPermut/lib/$(TARGET)
	rm -rf $(INSTALL_PATH)/hpcPermut/include/cuda/$(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET) *.optrpt
