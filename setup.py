# Suggestion for setup.py from
#   https://uoftcoders.github.io/studyGroup/lessons/python/packages/lesson/
# and
#   https://packaging.python.org/tutorials/packaging-projects/

import setuptools
import glob, re

with open("README.md", "r") as fh:
    long_description = fh.read()

def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
                       open(project + '/__init__.py').read())
    return result.group(1)


def get_version():
    """Get the version number of pyminbar"""
    # ## the original inspired by Paul's Aegean package
    # ## does not work with siplified imports trick in __init__.py
    import minbar
    return minbar.__version__

package_name='minbar'

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
    # see https://semver.org for guidelines
    # *strongly* suggested for sharing
    version=get_property('__version__', package_name),
    # The license can be anything you like
    license='GNU GPLv3',
    description='Python package for analysing observational data of thermonuclear (type-I) X-ray bursts observed by RXTE/PCA, BeppoSAX/WFC and INTEGRAL/JEM-X', 
    long_description=long_description,
    long_description_content_type="text/markdown",
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)

