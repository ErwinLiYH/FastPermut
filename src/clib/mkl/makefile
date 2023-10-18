#module load intel-20.4/cc
#module load intel-20.4/cmkl
#source python env, example: source /work/m22oc/m22oc/s2303098/pyenv4fyp/bin/activate

# compiler
CC := icpc
# compilation flags
CFLAGS := -Wall -Wextra -shared -fPIC -O3 -mkl -lmpi -qopenmp
# the target shared library
TARGET := libpermut.so
# source files
SRCS := permutation.cpp utils_mkl.cpp permutation_omp.cpp ../common/utils.cpp
# object files
OBJS := $(SRCS:.cpp=.o)
# install paths
INSTALL_PATH := $(PERMUT_PATH)
# header files
HEADERS := permutation.hpp

.PHONY: all install uninstall clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

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
