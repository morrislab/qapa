import os.path
import sys
from setuptools import find_packages
from distutils.core import setup

if sys.version_info[:2] < (3, 8):
    sys.stderr.write("Error: at least Python 3.8 or later is required.\n")
    sys.exit(1)

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'qapa/version.py')).read())

setup(name='qapa',
      version=__version__,  # noqa: F821
      description='RNA-seq Quantification of Alternative Polyadenylation (QAPA)',
      url='http://github.com/morrislab/qapa',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='GPLv3',
      packages=find_packages(),
      scripts=['scripts/create_merged_data.R',
               'scripts/compute_pau.R'],
      install_requires=['setuptools',
                        'pandas >= 0.24',
                        'numpy >= 1.10.0',
                        'biopython >= 1.76',
                        'pybedtools >= 0.7.9'
                        ],
      entry_points={
          'console_scripts': [
              'qapa = qapa.qapa:main'
          ]
      },
      python_requires='~=3.8',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8',
      ],
      zip_safe=False
      )
