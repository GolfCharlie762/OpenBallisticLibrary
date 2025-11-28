#!/usr/bin/env python3
"""
Setup file for Open Ballistic Library
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="open-ballistic-library",
    version="0.1.0",
    author="Open Ballistic Library Contributors",
    author_email="info@example.com",
    description="An open-source library for ballistic calculations and orbital mechanics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/example/open-ballistic-library",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Utilities",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
        "visualization": [
            "matplotlib>=3.0",
            "plotly>=4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "open-ballistic=open_ballistic_library.main:main",
        ],
    },
)