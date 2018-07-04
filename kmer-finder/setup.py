from setuptools import setup, find_packages
import os
import re

version = re.search(
    '__version__ = "([^\']+)"',
    open(
        os.path.join(os.path.dirname(__file__), 'kmerFinder/__init__.py')
    ).read()).group(1)

setup(
    name="KmerFinder",
    version=version,
    author="Your Name Here",
    author_email="youremail@example.com",
    url="https://github.com/yourusername/package-name-as-shown-on-pypi",
    description="Short package description",
    long_description='\n\n'.join((
        open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    )),
    entry_points={
        'console_scripts': [
            'findTemplate = kmerFinder.template.find:findTemplate',
            'maketemplatedb = kmerFinder.template.make:makeTemplateDB',
            'getTax = kmerFinder.output.taxonomy:getTaxonomy',
            'createTable = kmerFinder.output.table:createTSV',
            'readPrintKmerDB = kmerFinder.output.kmer:readPrintKmerDB',
            'makeTree = kmerFinder.output.tree:makeTree',
            'makeorganismDB = kmerFinder.template.organism:makeorganismDB',
        ]
    },
    packages=find_packages(),
    include_package_data=True,
    # test_suite = "kmerFinder.tests.make"
)
