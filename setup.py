# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


starcode_module=Extension('dnastorage.starcode_bindings',
                          sources = ['dnastorage/util/starcode_bindings.c', 'other_software/starcode/src/trie.c',
                                     'other_software/starcode/src/starcode.c'],
                          include_dirs = ['other_software/starcode/src'],
                          extra_compile_args = ['-std=c11',"-g3", "-O0", "-g3", "-O0" ]
                          )


fasthedges = Extension('dnastorage.codec.fasthedges',
                       sources = ['dnastorage/codec/fasthedges/module.cpp', \
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp'],
                      #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                      extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                       language='c++',)

codewordhedges = Extension('dnastorage.codec.codewordhedges',
                       sources = ['dnastorage/codec/fasthedges/codewordhedges_module.cpp',
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp',
                                  'dnastorage/codec/fasthedges/codeword_hedges.cpp'],
                          #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                           extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                           undef_macros=['NDEBUG'],
                           language='c++',)
hedges_hooks = Extension('dnastorage.codec.hedges_hooks',
                         sources = ['dnastorage/codec/fasthedges/hedges_hooks_module.cpp',
                                  'dnastorage/codec/fasthedges/fast_hedges.cpp'],
                         #extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O0",'-g3','-D DEBUG'],
                         extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],
                         undef_macros=['NDEBUG'],
                         language='c++',)

generate = Extension('dnastorage.util.generate',                                                                                                                           
                     sources = ['dnastorage/util/random_int.cpp'],                                                                                                        
                     extra_compile_args=["-std=c++11", "-Wall", "-Wextra","-O3"],                                                                                                  
                     undef_macros=['NDEBUG'],                                                                                                                                      
                     language='c++',)                                                                                                                                              

 
setup(
    name='dnastorage',
    version='1.1.1',
    description='DNA-based data storage modeling and simulation package',
    long_description=readme,
    author='James Tuck',
    author_email='jtuck@ncsu.edu',
    url='https://github.ncsu.edu/jtuck/',
    license=license,
    packages=find_packages(exclude=( 'tests','docs', 'tools', 'other_software')),
    ext_modules = [fasthedges,codewordhedges,generate,hedges_hooks,starcode_module]
)
