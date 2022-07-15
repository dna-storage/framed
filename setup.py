# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


starcode_module=Extension('starcode_bindings',
                          sources = ['dnastorage/util/starcode_bindings.c', 'other_software/starcode/src/trie.c',
                                     'other_software/starcode/src/starcode.c'],
                          include_dirs = ['other_software/starcode/src'],
                          extra_compile_args = ['-std=c99',"-g3", "-O0", "-g3", "-O0" ]
                          )


fasthedges = Extension('dnastorage.codec.fasthedges',
                       sources = ['dnastorage/codec/fasthedges/module.cpp', \
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp'],
                       extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                       #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                       language='c++',)

codewordhedges = Extension('dnastorage.codec.codewordhedges',
                       sources = ['dnastorage/codec/fasthedges/codewordhedges_module.cpp',
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp'],
                          extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                           #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                           undef_macros=['NDEBUG'],
                           language='c++',)




setup(name='generate', version = '1.0', ext_modules=[Extension('generate',['dnastorage/util/random_int.c'])])
setup(name='starcode_bindings', version='1.0', description='starcode bindings extensions',
      ext_modules=[starcode_module])

setup(
    name='dnastorage',
    version='1.1.0',
    description='DNA-based data storage modeling and simulation package',
    long_description=readme,
    author='James Tuck',
    author_email='jtuck@ncsu.edu',
    url='https://github.ncsu.edu/jtuck/',
    license=license,
    packages=find_packages(exclude=( 'tests','docs', 'tools', 'other_software')),
    ext_modules = [fasthedges,codewordhedges]
)
