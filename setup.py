from setuptools import setup, find_packages

setup(name="linnea",
    version="0.1",
    description="An experimental tool for the automatic generation of optimized code for linear algebra problems.",
    url="https://github.com/HPAC/linnea",
    author="Henrik Barthels",
    author_email="barthels@aices.rwth-aachen.de",
    license="GPLv3",
    packages=find_packages(exclude=('tests', )),
    zip_safe=True,
    python_requires=">=3.6",
    install_requires=[
        "matchpy >= 0.3",
    ],
    entry_points = {
        "console_scripts": ["linnea=linnea.__main__:main"],
    }
)
