from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import sys
import os

compile_cpu_in = "cpu" in sys.argv
compile_gpu_in = "gpu" in sys.argv

if os.environ.get("INSTALL_SELECTOR")==None or os.environ.get("INSTALL_SELECTOR")=="all":
    compile_cpu = True
    compile_gpu = True
    print("Install both CPU and GPU version")
elif os.environ.get("INSTALL_SELECTOR")=="cpu":
    compile_cpu = True
    compile_gpu = False
    print("Only install CPU version")
elif os.environ.get("INSTALL_SELECTOR")=="gpu":
    compile_cpu = False
    compile_gpu = True
    print("Only install GPU version")
else:
    raise Exception("INSTALL_SELECTOR should be in [\"all\", \"cpu\", \"gpu\"], or do not define it for \"all\". Got %s."%os.environ.get("INSTALL_SELECTOR"))

if compile_cpu_in:
    sys.argv.remove("cpu")
if compile_gpu_in:
    sys.argv.remove("gpu")

extensions = []

if compile_cpu:
    cpu_extension = Extension(
                        "fastpermut.permutation",  # Name of the resulting python extension
                        sources=["src/sapermut/pypermutation.pyx"],  # Cython source file
                        libraries=["permut"],  # Name of the C library without the 'lib' prefix and '.so' suffix
                        library_dirs=[os.path.expanduser("~/hpcPermut/lib"),
                                    os.path.join(os.environ.get('PERMUT_PATH'),"hpcPermut/lib")],  # Path to the directory of your .so file
                        language='c++',
                        include_dirs=[numpy.get_include(),
                                    os.path.expanduser("~/hpcPermut/include"),
                                    os.path.join(os.environ.get('PERMUT_PATH'),"hpcPermut/include/mkl")],
                        # Path to the directory of your header files
                )
    extensions.append(cpu_extension)
if compile_gpu:
    gpu_extension = Extension(
                        "fastpermut.permutation_cuda",  # Name of the resulting python extension
                        sources=["src/sapermut/pypermutation_cuda.pyx"],  # Cython source file
                        libraries=["permut_cuda"],  # Name of the C library without the 'lib' prefix and '.so' suffix
                        library_dirs=[os.path.expanduser("~/hpcPermut/lib"),
                                    os.path.join(os.environ.get('PERMUT_PATH'),"hpcPermut/lib")],  # Path to the directory of your .so file
                        language='c++',
                        include_dirs=[numpy.get_include(),
                                    os.path.expanduser("~/hpcPermut/include"),
                                    os.path.join(os.environ.get('PERMUT_PATH'),"hpcPermut/include/cuda"),
                                    "/mnt/lustre/indy2lfs/sw/nvidia/hpcsdk-2211/Linux_x86_64/22.11/cuda/11.8/targets/x86_64-linux/include"], #used in cirrus, to find driver_types.h
                        # Path to the directory of your header files
                )
    extensions.append(gpu_extension)

setup(
    name="Fast permutation",
    version="0.0.5",
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    author="Yuanhao Li",
    author_email="erwinli@qq.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires="==3.9.*",
    ext_modules=cythonize(extensions, language_level=3),
)
