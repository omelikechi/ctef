from os import path
from setuptools import find_packages, setup

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="nflows",
    version="0.1dev",
    description="Cayley transform ellipsoid fitting",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/omelikechi/ctef",
    author="Omar Melikechi",
    packages=find_packages(exclude=["tests"]),
    license="MIT",
    install_requires=[
        "matplotlib",
        "numpy",
        "tqdm",
    ],
    dependency_links=[],
)