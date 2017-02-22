from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('version.txt').read()) # loads __version__

setup(name='refineomatic',
      version=__version__,
      author='Joel Boyd',
      description='Simple pipeline for lazy people who want refine their bins',
      long_description=readme,
      license='GPL3+',
      keywords="",
      packages=find_packages(exclude='docs'),
      install_requires=('biopython >=1.64'),
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      url='https://github.com/geronimp/refineomatic',
      scripts=['bin/refine-o-matic.py'],
      data_files=[])
