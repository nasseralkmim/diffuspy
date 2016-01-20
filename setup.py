from setuptools import setup

setup(name='poisson',
      version='0.1',
      description='Poisson equation solver',
      url='nasseralkmim.github.io',
      author='Nasser Alkmim',
      author_email='nasser.alkmim@gmail.com',
      license='MIT',
      packages=['poisson'],
      install_requires=[
          'numpy',
          'matplotlib',
          'networkx',
          'scipy'
      ],
      zip_safe=False)
