# minimal setup.py to be able to use the -e flag (pip install -e .)

from setuptools import find_packages, setup

setup(
    package_data={"pyranges": ["data/*.bam", "data/ncbi.*", "data/bigwig.bw"]},
    include_package_data=True,
    packages=find_packages(),
)
