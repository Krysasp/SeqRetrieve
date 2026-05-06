#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name='SeqRetrieve',
    version='1.0.0',
    description='GenBank sequence retrieval and analysis toolkit',
    author='Your Name',
    author_email='jonatasp92@gmail.com',
    url='https://github.com/Krysasp/RetrieveSeq',
    packages=find_packages(),
    scripts=['bin/retrieveseq'],
    install_requires=[
        'biopython>=1.80',
        'numpy>=1.20',
    ],
    python_requires='>=3.7',
)
