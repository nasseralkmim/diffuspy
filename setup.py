from setuptools import setup

setup(name='diffuspy',
      version='0.1',
      description='General diffusion equation in 2d solver using FEM',
      url='nasseralkmim.github.io',
      author='Nasser Alkmim',
      author_email='nasser.alkmim@gmail.com',
      license='MIT',
      packages=['diffuspy'],
      install_requires=[
          'numpy',
          'matplotlib',
          'networkx',
          'scipy'
      ],
      zip_safe=False)
