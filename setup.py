import os
from setuptools import setup


def read(fname):
    """
    Utility function to read information from files.
    Used to read README for the long_description and VERSION for the version number.
    It's nice, because now 1) we have a top level
    README file 2) it's easier to type in the README file than to put a raw
    string in below and 3) we can use a central place to keep version information.
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="pyrates",
    version=read('VERSION'),
    author="Peter Humburg, John Judge, David Kohn and Team Pirates",
    author_email="p.humburg@garvan.org.au",
    description=("Reduces errors in deep sequencing data by computing consensus",
                 "sequences for duplicate reads"),
    license="MIT",
    keywords="sequencing error-correction consensus fastq",
    url="https://github.com/humburg/pirates",
    packages=['pyrates'],
    entry_points={
        "console_scripts":['pyrates = pyrates.cmd_consensus:main']
    },
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
    ],
)
