from skbuild import setup
from setuptools import find_packages

try:
    import numpy

    def get_include(*args, **kwargs):
        return numpy.get_include()


except ImportError:

    def get_include(*args, **kwargs):
        import numpy

        return numpy.get_include()



if __name__ == "__main__":
    setup(
        name="libvlmc",
        version="0.1",
        description="Genomic signatures and related distance measures for genomic sequences.",
        packages=find_packages(),
        install_requires=['cython', 'numpy'],
        include_dirs=[get_include()],
        setup_requires=['numpy'],
        zip_safe=False,
        classifier=["Private :: Do Not Upload"],
        cmake_args=['-DCMAKE_BUILD_TYPE=Release'],
        python_requires='>=3',
    )
