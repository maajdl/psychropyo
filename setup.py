# from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name='psychropyo',
    version='0.0.1',
    description='psychrometry using the pyomo solver',
    author='maajdl',
    author_email='',
    packages=find_packages(),
    install_requires=['pyomo','holoviews','pandas','numpy'],
    license='',
)
