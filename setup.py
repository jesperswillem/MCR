import re
from setuptools import setup, find_packages

version = "0.2.0"

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

setup(
    name='uppMCR',
    version='0.2',
    zip_safe=True,
    packages=['MCR',],
    scripts=["./bin/run_mcr.py", "./bin/mcr_plot.py", "./bin/construct_clusters.py"],
    test_suite="tests",
    url='https://github.com/jesperswillem/MCR',
    license='All rights reserved. (for now)',
    author='Florian van der Ent, Willem Jespers',
    author_email='florianvdent@gmail.com',
    long_description=long_descr,
    description='Set of tools and a library for the generation of virtual screening libraries from multi-component reactions.'
)
