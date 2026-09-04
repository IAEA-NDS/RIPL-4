
# OS
import os

# Setuptools
import setuptools


def read(rel_path: str) -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, rel_path), 'r', encoding='utf-8') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


# Grab contents of readme Markdown file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Requirements (one per line; skip blanks and comments)
with open("requirements.txt", "r") as fh:
    reqs = [
        line.strip()
        for line in fh
        if line.strip() and not line.strip().startswith("#")
    ]

# Optional features. pandas/numpy are only needed for their respective
# exporters (Database.to_dataframe / to_numpy and their riplpy.* wrappers);
# the package imports and the CSV/JSON/list exporters work without either.
extras = {
    "dataframe": ["pandas"],
    "numpy": ["numpy"],
}
extras["all"] = sorted({dep for deps in extras.values() for dep in deps})


setuptools.setup(
    name="Reference-Input-Parameter-Library",
    version=get_version("riplpy/__init__.py"),
    author="Matthew Mumpower",
    author_email="matthew@mumpower.net",
    description="A Python package for the Reference Input Parameter Library (RIPL)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mrmumpow/riplpy",
    project_urls={
        "Source": "https://github.com/mrmumpow/riplpy",
        "Bug Tracker": "https://github.com/mrmumpow/riplpy/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    # Explicit package discovery (restricted to the riplpy tree) so installs
    # are deterministic and never sweep in tests/ under flat-layout
    # auto-discovery.
    packages=setuptools.find_packages(include=["riplpy", "riplpy.*"]),
    include_package_data=True,
    python_requires=">=3.10",
    install_requires=reqs,
    extras_require=extras,
)
