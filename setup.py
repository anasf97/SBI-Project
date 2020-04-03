from setuptools import setup, find_packages
import os

with open("README.md", encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'ComplexAssembler',
    version = '1.0',
    description = "Macrocomplex builder",
    scripts = ['ComplexAssembler/assemble.py'],
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    entry_points={
        'console_scripts': [
            'assemble = ComplexAssembler.assemble:main'
        ]
    },
    url = 'https://github.com/anasf97/SBI-Project',
    author = 'Eva Brigos Barril, Jordi Busoms Pratdesaba, Ana Sanchez Fernandez',
    classifiers =[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License"
    ],
    packages = find_packages(),
    install_requires = ['biopython', 'matplotlib', 'numpy', 'scipy']
)
