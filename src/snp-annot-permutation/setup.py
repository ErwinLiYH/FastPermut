from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os

clib_path="" # c library path

extensions = [
    Extension(
        "sapermut.permutation",  # Name of the resulting python extension
        sources=["src/sapermut/pypermutation.pyx"],  # Cython source file
        libraries=["permut"],  # Name of the C library without the 'lib' prefix and '.so' suffix
        library_dirs=[os.path.expanduser("~/hpcPermut/lib")],  # Path to the directory of your .so file
        language='c++',
        include_dirs=[numpy.get_include(),
                      os.path.expanduser("~/hpcPermut/include")],
        # Path to the directory of your header files
    ),
    Extension(
        "sapermut.permutation_cuda",  # Name of the resulting python extension
        sources=["src/sapermut/pypermutation_cuda.pyx"],  # Cython source file
        libraries=["permut_cuda"],  # Name of the C library without the 'lib' prefix and '.so' suffix
        library_dirs=[os.path.expanduser("~/hpcPermut/lib")],  # Path to the directory of your .so file
        language='c++',
        include_dirs=[numpy.get_include(),
                      os.path.expanduser("~/hpcPermut/include"),
                      "/mnt/lustre/indy2lfs/sw/nvidia/hpcsdk-2211/Linux_x86_64/22.11/cuda/11.8/targets/x86_64-linux/include"],
        # Path to the directory of your header files
    )
]
setup(
    name="SNP-annot-permutation",
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
