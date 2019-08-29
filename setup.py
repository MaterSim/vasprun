from distutils.core import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="vasprun",
    version="0.0dev",
    author="Qiang Zhu",
    author_email="qiang.zhu@unlv.edu",
    description="Python code for vasprun.xml analyais",
    long_description=long_description,
    url="https://github.com/qzhu2017/vasprun",
    packages=['vasprun'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    requires=['pymatgen', 'numpy', 'scipy', 'matplotlib'],
)
