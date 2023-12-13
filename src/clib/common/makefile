CC := g++
# compilation flags
CFLAGS := -Wall -Wextra -shared -fPIC -O3
# the target shared library
TARGET := libpermut_utils.so
# source files
SRCS := utils.cpp
# object files
OBJS := $(SRCS:.cpp=.o)
# install paths
ifndef PERMUT_PATH
	INSTALL_PATH := $(HOME)
else
	INSTALL_PATH := $(PERMUT_PATH)
endif
# header files
HEADERS := utils.hpp config.hpp

.PHONY: all install uninstall clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

install:
	mkdir -p $(INSTALL_PATH)/hpcPermut/lib
	cp $(TARGET) $(INSTALL_PATH)/hpcPermut/lib

	mkdir -p $(INSTALL_PATH)/hpcPermut/include/common
	cp $(HEADERS) $(INSTALL_PATH)/hpcPermut/include/common

uninstall:
	rm -rf $(INSTALL_PATH)/hpcPermut/lib/$(TARGET)
	$(foreach header,$(HEADERS),rm -f $(INSTALL_PATH)/hpcPermut/include/common/$(header);)

clean:
	rm -f $(OBJS) $(TARGET) *.optrpt
