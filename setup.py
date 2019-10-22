from distutils.core import setup
import setuptools  # noqa

setup(
    name="vasprun-xml",
    version="1.0",
    author="Qiang Zhu",
    author_email="qiang.zhu@unlv.edu",
    description="description",
    url="https://github.com/qzhu2017/vasprun/archive/v1.0.zip",
    packages=['vasprun'],
    scripts=['vasprun/vasprun'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'pymatgen>=2018.12.12', 
        'lxml>=4.2.5', 
        'pandas>=0.23.4', 
        'numpy>=1.13.3', 
        'scipy>=1.1.0', 
        'matplotlib>=2.0.0'],
    python_requires='>=3.6.1',
    license='MIT',
)
