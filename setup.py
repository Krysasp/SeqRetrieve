#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name='SeqRetrieve',
    version='1.0.0',
    description='Comprehensive toolkit for retrieving, sorting, aligning, filtering, and reporting GenBank sequences',
    author='RetrieveSeq Contributors',
    author_email='jonatasp92@gmail.com',
    url='https://github.com/Krysasp/SeqRetrieve',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    scripts=['bin/retrieveseq'],
    install_requires=[
        'biopython>=1.80',
        'numpy>=1.20',
    ],
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Science/Research :: Bioinformatics',
    ],
    keywords='genbank sequence retrieval metadata alignment filtering bioinformatics',
)
