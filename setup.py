from setuptools import setup
from distutils.extension import Extension
import numpy
from Cython.Build    import cythonize
from Cython.Compiler import Options as CCoptions

CCoptions.annotate = True

extension = Extension("transformations._ctransformations",
                      sources=["transformations/_ctransformations.pyx"],
                      include_dirs=[numpy.get_include()],
                      extra_compile_args=["-O3"],
                      language="c++",
                      define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name='transformations',
    version='0.0.1',
    description='Trajectory transformations for MDAnalysis',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    packages=["transformations"],
    author='Santeri Paajanen',
    author_email='santeri.e.paajanen@helsinki.fi',
    ext_modules = cythonize(extension, language_level = "3")
)