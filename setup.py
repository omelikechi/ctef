from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="ctep",
    version="0.0.1",
    description="Cayley transform ellipsoid fitting",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/omelikechi/ctef",
    author="Omar Melikechi",
    license="MIT",
    install_requires=[
        "matplotlib",
        "numpy",
        "tqdm",
    ],
    dependency_links=[],
)