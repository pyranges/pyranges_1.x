# minimal setup.py to be able to use the -e flag (pip install -e .)

from setuptools import setup, find_packages

setup(
    package_data={"pyranges": ["data/*.bam"]},
    include_package_data=True,
    packages=find_packages(),
)
