# -*- coding: utf-8 -*-

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


fasthedges = Extension('dnastorage.codec.fasthedges',
                       sources = ['dnastorage/codec/fasthedges/module.cpp', \
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp'],
                      #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                      extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                       language='c++',)

generate = Extension('dnastorage.util.generate',                                                                                                                           
                     sources = ['dnastorage/util/random_int.cpp'],                                                                                                        
                     extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],                                                                                                  
                     undef_macros=['NDEBUG'],                                                                                                                                      
                     language='c++',)                                                                                                                                              


class OptionalBuildExt(build_ext):
    """Build C++ extensions when possible; skip gracefully if compilation fails.

    In environments without a suitable C++ compiler the pure-Python fallbacks
    (e.g. ``dnastorage/util/generate.py``) will be used automatically.
    """

    def run(self):
        try:
            super().run()
        except Exception as exc:
            print(
                f"WARNING: C++ extension build step failed ({exc}). "
                "The package will still work using pure-Python fallbacks, "
                "but performance-sensitive features (fasthedges, generate) "
                "will be slower."
            )

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as exc:
            print(
                f"WARNING: Could not build extension '{ext.name}': {exc}. "
                "Falling back to the pure-Python implementation."
            )


setup(
    name='dnastorage',
    version='1.1.1',
    description='DNA-based data storage modeling and simulation package',
    long_description=readme,
    author='James Tuck, Kevin Volkel',
    author_email='jtuck@ncsu.edu, kvolkel@ncsu.edu',
    url='',
    license=license,
    packages=find_packages(exclude=( 'tests','docs', 'tools', 'other_software')),
    ext_modules = [fasthedges,generate],
    cmdclass={'build_ext': OptionalBuildExt},
)
