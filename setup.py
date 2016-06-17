from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import cython_gsl
import numpy as np
from distutils.core import setup

extensions = [
    Extension("flow_pop.cy_particles",
              sources=["flow_pop/cy_particles.pyx"],
              language="c", libraries = cython_gsl.get_libraries(),
              library_dirs = [cython_gsl.get_library_dir()],
              include_dirs = [cython_gsl.get_cython_include_dir(), np.get_include()])
]


setup(
    name='particle_populations_in_flow',
    version='0.1',
    packages=['flow_pop'],
    include_dirs = [cython_gsl.get_include(), np.get_include()],
    ext_modules = cythonize(extensions, annotate=True),
    url='',
    license='',
    author='bryan',
    author_email='bweinstein@seas.harvard.edu',
    description='',
)
