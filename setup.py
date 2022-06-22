from setuptools import setup, Extension
import numpy
from Cython.Build    import cythonize
#from Cython.Compiler import Options as CCoptions

#CCoptions.annotate = True

extension = Extension("transformations._ctransformations",
                      sources=["lib/_ctransformations.pyx"],
                      include_dirs=[numpy.get_include()],
                      extra_compile_args=["-O3"],
                      language="c++",
                      define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])


setup(
    packages=["transformations"],
    ext_modules = cythonize(extension, language_level = "3")
)