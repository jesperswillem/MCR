import re
from setuptools import setup, find_packages

version = "0.3.0"

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

setup(
    name='MCR_devel',
    version='0.3',
    zip_safe=True,
    packages=['MCR',],
    include_package_data=True,
    scripts=["./bin/mcr_run.py", "./bin/mcr_plot.py", "./bin/mcr_reconstruct_clusters.py", "./bin/mcr_find_hits.py"],
    test_suite="tests",
    url='https://github.com/jesperswillem/MCR',
    license='All rights reserved. (for now)',
    author='Willem Jespers, Florian van der Ent, Hugo Gutierrez de Teran',
    author_email='florianvdent@gmail.com',
    long_description=long_descr,
    description='Set of tools and a library for the generation of virtual screening libraries from multi-component reactions.'
)
