# Suggestion for setup.py from
#   https://uoftcoders.github.io/studyGroup/lessons/python/packages/lesson/
# and
#   https://packaging.python.org/tutorials/packaging-projects/

import setuptools
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='minbar',
    url='https://github.com/outs1der/minbar',
    author='Duncan Galloway & Laurens Keek',
    author_email='Duncan.Galloway@monash.edu',
    # Needed to actually package something
    packages=['minbar'],
    # Needed for dependencies; probably not complete
    install_requires=['numpy', 'pandas', 'astropy', 'scipy'],
    # *strongly* suggested for sharing
    version='1.0.0',
    # The license can be anything you like
    license='GNU GPLv3',
    description='Python package for analysing observational data of thermonuclear (type-I) X-ray bursts observed by RXTE/PCA, BeppoSAX/WFC and INTEGRAL/JEM-X', 
    long_description=long_description,
    long_description_content_type="text/markdown",
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)

