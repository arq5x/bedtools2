import os
from setuptools import setup, find_packages

setup(
    name='bedtools',
    version="2.16.2",
    description='A flexible toolset for genome arithmetic',
    author='Aaron Quinkan',
    author_email='arq5x@virginia.edu',
    url='http://bedtools.googelcode.come',
    packages=find_packages(),
    include_package_data=True,
    #install_requires=libraries,
    zip_safe=False,
    #test_suite='readthedocs.tests.rtd_tests.runtests'
)