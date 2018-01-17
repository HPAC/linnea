from setuptools import setup

setup(name="linnea",
    version="0.1",
    # description="The funniest joke in the world",
    # url="http://github.com/storborg/funniest",
    author="Henrik Barthels",
    author_email="barthels@aices.rwth-aachen.de",
    # license="MIT",
    packages=["linnea"],
    zip_safe=True,
    python_requires=">=3.5.3",
    install_requires=[
        "matchpy >= 0.3",
    ],
    entry_points = {
        "console_scripts": ["linnea=linnea.__main__:main"],
    }
)