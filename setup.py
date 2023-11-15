# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension


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
    ext_modules = [fasthedges,generate]
)
