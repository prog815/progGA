from setuptools import setup, find_packages
from os.path import join, dirname

import progGA

setup(
    name='progGA',
    version='0.0.7',
    description='Библиотека генетических алгоритмов',
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.md')).read(),
    include_package_data=True,
    author_email='eavprog@gmail.com'
)
