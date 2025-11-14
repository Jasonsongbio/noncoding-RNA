from setuptools import find_packages, setup

setup(
    name="phylo-pruner",
    version="0.1.0",
    description="Annotate and prune phylogenetic trees using the NCBI taxonomy database.",
    author="Your Name",
    license="MIT",
    packages=find_packages(include=["phylo_pruner", "phylo_pruner.*"]),
    install_requires=[
        "ete3",
    ],
    entry_points={
        "console_scripts": [
            "phylo_pruner=phylo_pruner.main:cli_entry",
        ]
    },
    python_requires=">=3.8",
    include_package_data=True,
)
