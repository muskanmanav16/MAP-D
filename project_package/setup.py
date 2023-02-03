#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = [
    'numpy',
    'Click',
    'networkx',
    'matplotlib',
    'pandas',
    'requests',
    'tqdm',
    'scipy',
    'sqlalchemy',
    'sqlalchemy_utils',
    'xmltodict',
    'fastapi',
    'fastapi-pagination',
    'python-multipart',
    'uvicorn',
    'flask',
    'PyMySQL',
    'cryptography',
    'scipy',
    'spacy>=3.4.0,<3.5.0',
    'scispacy==0.5.1',
    'en-ner-bionlp13cg-md @ https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_bionlp13cg_md-0.5.1.tar.gz',
    # 'spacy',
    'tensorflow',



]

test_requirements = ['pytest>=3', ]

setup(
    author="Muskan Astha Deepika Parinishtha",
    author_email='mapd@gmx.net',
    python_requires='>=2.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Map-D Package",
    entry_points={
        'console_scripts': [
            'mapd=mapd.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='mapd',
    name='mapd',
    packages=find_packages(include=['mapd', 'mapd.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
