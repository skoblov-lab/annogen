from distutils.core import setup
import os

from Cython.Build import cythonize

os.environ['CFLAGS'] = '-O3 -Wall -std=c++11 -stdlib=libc++'

setup(
    name="annogen",
    ext_modules=cythonize("annogen/*.pyx",
                          include_path=["annogen/"],
                          language="c++"),
    packages=["annogen"],
    install_requires=["cython>=0.27"]
)
