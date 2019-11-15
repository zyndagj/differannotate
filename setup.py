try:
	from setuptools import setup
except:
	from distutils.core import setup
from glob import glob

setup(name = "differannotate",
	version = "0.0.1",
	author = "Greg Zynda",
	author_email="gzynda@tacc.utexas.edu",
	license="BSD-3",
	description="A tool for comparing GFF3 annotations",
	install_requires=["quicksect","numpy","matplotlib-venn"],
	tests_requires=["quicksect","numpy","matplotlib-venn"],
	packages = ["differannotate"],
	entry_points = {'console_scripts': ['differannotate=differannotate:main']},
	options = {'build_scripts': {'executable': '/usr/bin/env python'}},
	test_suite = "tests")
	#scripts = list(glob('scripts/*')),
