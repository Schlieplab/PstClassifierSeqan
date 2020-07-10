from skbuild import setup
from setuptools import find_packages

if __name__ == "__main__":
    setup(
        name="libvlmc",
        version="0.1",
        description="Genomic signatures and related distance measures for genomic sequences.",
        packages=find_packages(),
        install_requires=['cython'],
        zip_safe=False,
        classifier=["Private :: Do Not Upload"],
        cmake_args=['-DCMAKE_BUILD_TYPE=Release'],
        python_requires='>=3',
    )
