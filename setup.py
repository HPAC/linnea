from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="linnea",
    version="0.1.dev1",
    description="An experimental tool for the automatic generation of optimized code for linear algebra problems.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HPAC/linnea",
    author="Henrik Barthels",
    author_email="barthels@aices.rwth-aachen.de",
    license="GPLv3",
    packages=find_packages(exclude=('tests', )),
    include_package_data=True,
    zip_safe=True,
    python_requires=">=3.6",
    install_requires=[
        "matchpy >= 0.5.2",
        "tatsu >= 4.4"
    ],
    entry_points = {
        "console_scripts": ["linnea=linnea.__main__:main"],
    },
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent"
    ]
)
