import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="agamPrimer",
    version="0.1.4",
    author="Sanjay Curtis Nagi",
    author_email="sanjay.c.nagi@gmail.com",
    description="A small package to store some functions for the PrimerDesignAg notebook, to improve readability.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sanjaynagi/PrimerDesignAg",
    project_urls={
        "Bug Tracker": "https://github.com/sanjaynagi/PrimerDesignAg/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
   # package_dir={"": "agamPrimer"},
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
)
