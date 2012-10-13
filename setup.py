from distutils.core import setup

from version import __version__

setup(name='hyperfreq',
        version=__version__,
        url='http://github.com/fhcrc/hyperfreq',
        description='Hypermutation evaluation software',
        author="Christopher Small",
        author_email="csmall@fhcrc.org",
        scripts=['hyperfreq/scripts/hyperfreq'],
        packages=['hyperfreq'],
        requires=['sekhon', 'biopython'])

