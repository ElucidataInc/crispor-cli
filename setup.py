
from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='crispor_cli',
    version='0.1.0',
    description='Command line tool for crispr offtarget finding extracted from crispor website',
    long_description=readme,
    author='Kailash Yadav',
    author_email='kailash.yadav@elucidata.io',
    url='https://github.com/ElucidataInc/crispor-cli',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)