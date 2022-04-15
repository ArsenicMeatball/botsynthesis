from setuptools import setup, find_packages

setup(name='botsynthesis',
      version='4.0.1',
      author='Jean-David Rousseau',
      author_email='arsenicmeatball@yahoo.com',
      packages=find_packages(exclude=['tests', 'tests.*']),
      url='https://github.com/ArsenicMeatball/BOTS_development',
      scripts=[],
      license='LICENSE.md',
      description='Protein Sequence Optimizer',
      long_description=open('README.md').read(),
      install_requires=[
          'python_codon_tables',
          'biopython'
      ],
      )
