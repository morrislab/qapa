import sys
from setuptools import find_packages
from distutils.core import setup

if sys.version_info[:2] < (2, 7):
    sys.stderr.write("At least Python 2.7, or Python 3.5 or later "
                     "is required\n")
    sys.exit(1)
elif sys.version_info[0] == 3 and sys.version_info[:2] < (3, 5):
    sys.stderr.write("At least Python 3.5 or later is required.\n")
    sys.exit(1)


setup(name='qapa',
      version='1.0.0',
      description='RNA-seq Quantification of Alternative Polyadenylation (QAPA)',
      url='http://github.com/morrislab/qapa',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='See LICENSE',
      packages=find_packages(),
      scripts=['scripts/create_merged_data.R',
               'scripts/compute_pau.R'],
      install_requires=['setuptools',
                        'pandas >= 0.17',
                        'numpy >= 1.10.0',
                        'biopython >= 1.66',
                        'pybedtools >= 0.7.9'],
      entry_points={
            'console_scripts': [
                'qapa = qapa.qapa:main'
            ]
      },
      zip_safe=False
      )
