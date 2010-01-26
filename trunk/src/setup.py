__author__="Halikov_AA"
__date__ ="$26.01.2010 14:27:30$"

from setuptools import setup,find_packages

setup (
  name = 'pystab',
  version = '0.1',
  packages = find_packages(),

  # Declare your packages' dependencies here, for eg:
  install_requires=['foo>=3'],

  # Fill in these to make your Egg ready for upload to
  # PyPI
  author = 'Artur',
  author_email = 'akhalikoff@gmail.com',

  summary = '',
  url = '',
  license = '',
  long_description= 'Long description of the package',

  # could also include long_description, download_url, classifiers, etc.

  
)