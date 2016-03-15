from setuptools import setup

setup(name='pyheat',
      version='0.1',
      description='heat equation solver',
      url='nasseralkmim.github.io',
      author='Nasser Alkmim',
      author_email='nasser.alkmim@gmail.com',
      license='MIT',
      packages=['pyheat'],
      install_requires=[
          'numpy',
          'matplotlib',
          'networkx',
          'scipy'
      ],
      zip_safe=False)
