from distutils.core import setup
import setuptools  # noqa
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

exec(open('vasprun/version.py').read())

setup(
    name="vasprun-xml",
    version=__version__,
    author="Qiang Zhu",
    author_email="qiang.zhu@unlv.edu",
    description="A python package for quick analysis of vasp calculation",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url="https://github.com/qzhu2017/vasprun",
    packages=['vasprun'],
    scripts=['scripts/vasprun'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'lxml>=4.2.5', 
        'pandas>=0.23.4', 
        'numpy>=1.13.3', 
        'scipy>=1.1.0', 
        'matplotlib>=2.0.0',
        'tabulate'],
    python_requires='>=3.6.1',
    license='MIT',
)
