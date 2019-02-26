from distutils.core import setup

from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [Extension("my_bliss", ["my_bliss.pyx"],
                         include_dirs="bliss",
                         libraries=["bliss"])
]

setup(
    ext_modules=cythonize(ext_modules)
)
