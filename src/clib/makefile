ifndef PERMUT_PATH
	INSTALL_PATH := $(HOME)
else
	INSTALL_PATH := $(PERMUT_PATH)
endif

.PHONY: cpu gpu all common mkl cuda install_cpu install_gpu install uninstall uninstall_common uninstall_mkl uninstall_cuda clean

all: common mkl cuda

cpu: common mkl

gpu: common cuda

install_cpu: install_common install_mkl

install_gpu: install_common install_cuda

install: install_common install_mkl install_cuda

uninstall_cpu: uninstall_common uninstall_mkl

uninstall_gpu: uninstall_common uninstall_cuda

uninstall:
	rm -rf $(INSTALL_PATH)/hpcPermut

common:
	$(MAKE) -C common

mkl:
	$(MAKE) -C mkl

cuda:
	$(MAKE) -C cuda

install_common:
	$(MAKE) -C common install

install_mkl:
	$(MAKE) -C mkl install

install_cuda:
	$(MAKE) -C cuda install

uninstall_common:
	$(MAKE) -C common uninstall

uninstall_mkl:
	$(MAKE) -C mkl uninstall

uninstall_cuda:
	$(MAKE) -C cuda uninstall

clean:
	$(MAKE) -C common clean
	$(MAKE) -C mkl clean
	$(MAKE) -C cuda clean