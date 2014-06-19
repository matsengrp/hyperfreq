from setuptools import setup, find_packages

from hyperfreq.version import __version__

setup(name='hyperfreq',
        version=__version__,
        url='http://github.com/fhcrc/hyperfreq',
        description='Hypermutation evaluation software',
        author="Christopher Small",
        author_email="csmall@fhcrc.org",
        entry_points={
            'console_scripts': [
                'hyperfreq = hyperfreq.scripts.cli:main'
            ]},
        packages=find_packages(exclude=['tests']),
        install_requires=['betarat >=1.0, <2.0', 'alnclst >=0.0.1', 'biopython', 'fisher'])

