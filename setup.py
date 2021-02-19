"""Setup script for bfx_tools package."""

import sys

from setuptools import setup, find_packages

if sys.version_info[:2] < (3, 8):
    sys.exit('Sorry, Python 3.8 or newer is required.')


install_requires = [
    'requests',
    'setuptools_scm'
]


tests_require = [
    'coverage',
    'flake8',
    'flake8-bugbear',
    'flake8-builtins',
    'flake8-comprehensions',
    'flake8-blind-except',
    'flake8-docstrings',
    'flake8-mutable',
    'flake8-rst-docstrings',
    'mypy',
    'pyflakes',
    'pytest',
    'pytest-cache',
    'pytest-cov',
    'pytest-flake8',
    'pytest-html',
    'pytest-mypy',
    'pytest-xdist',
    'requests',
]


if __name__ == '__main__':
    setup(
        name='bfx_tools-Larry_Clos',
        version='0.1.0',
        description='Some Bioinformatics tools I developed in my own time.',
        author='Larry Clos',
        author_email='drlclos@gmail.com',
        license='Open',
        classifiers=[
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: OS Independent'
        ],
        packages=find_packages(),
        install_requires=install_requires,
        tests_require=tests_require,
        entry_points={},
        python_requires='>=3.8'
    )
