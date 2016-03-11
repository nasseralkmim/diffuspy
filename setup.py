from setuptools import setup

setup(name='pyisson',
      version='0.1',
      description='Poisson equation solver',
      url='nasseralkmim.github.io',
      author='Nasser Alkmim',
      author_email='nasser.alkmim@gmail.com',
      license='MIT',
      packages=['pyisson'],
      install_requires=[
          'numpy',
          'matplotlib',
          'networkx',
          'scipy'
      ],
      zip_safe=False)
