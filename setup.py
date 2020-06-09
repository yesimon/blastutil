# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='blastutil',
    version='0.1.0',
    description='Toolkit and library for working with BLAST',
    long_description=readme,
    author='Simon Ye',
    author_email='mail@yesimon.com',
    url='https://github.com/yesimon/blastutil',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    entry_points = {
        'console_scripts': [
            'blastutil = blastutil.__main__:main'
        ]
    }
)
