
import sys

from setuptools import setup, Extension, find_packages

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    raise SystemExit("PhasingTools requires Python 2.7")

__VERSION__ = "0.1.0"

DESC = 'Tools for analyzing complex HLA data from SMRT Sequencing'

scripts = [
    "bin/LociAnalysis",
]

required = [
    "numpy >= 1.10.1",
    "pbcore >= 1.2.6",
]

setup(
    name = 'LociAnalysis',
    version=__VERSION__,
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    url='https://github.com/bnbowman/LociAnalysis',
    description=DESC,
    license=open('LICENSES.txt').read(),
    packages = find_packages(),
    #package_dir = {'':'LociAnalysis'},
    zip_safe = False,
    scripts = scripts,
    install_requires = required,
)
