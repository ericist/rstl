import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rstl",
    version="0.1.3",
    author="Eric Rauch",
    author_email="ericist@pm.me",
    description="A Python port of R's stl function",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ericist/rstl",
    packages=["rstl"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy"]
)
