#!/usr/bin/env python
# encoding: utf-8

import re
import os

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext as _build_ext
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command.build_ext import build_ext as _build_ext


def find_boost(hint=None, verbose=True):
    """
    Find the location of the Boost include directory. This will return
    ``None`` on failure.

    """
    # List the standard locations including a user supplied hint.
    search_dirs = [] if hint is None else hint
    search_dirs += [
        "/usr/local/include",
        "/usr/local/homebrew/include",
        "/opt/local/var/macports/software",
        "/opt/local/include",
        "/usr/include",
        "/usr/include/local",
    ]

    # Loop over search paths and check for the existence of the required
    # header.
    for d in search_dirs:
        path = os.path.join(d, "boost", "math", "special_functions",
                            "ellint_3.hpp")
        if os.path.exists(path):
            # Determine the version.
            vf = os.path.join(d, "boost", "version.hpp")
            if not os.path.exists(vf):
                continue
            src = open(vf, "r").read()
            v = re.findall("#define BOOST_LIB_VERSION \"(.+)\"", src)
            if not len(v):
                continue
            v = v[0]
            if verbose:
                print("Found Boost version {0} in: {1}".format(v, d))
            return d
    return None


class build_ext(_build_ext):
    """
    A custom extension builder that finds the include directories for Boost.

    """

    def build_extension(self, ext):
        dirs = ext.include_dirs + self.compiler.include_dirs

        # Look for the Boost headers and make sure that we can find them.
        boost_include = find_boost(hint=dirs)
        if boost_include is None:
            raise RuntimeError("Required library Boost not found. "
                               "Check the documentation for solutions.")

        # Update the extension's include directories.
        ext.include_dirs += [boost_include]

        # Run the standard build procedure.
        _build_ext.build_extension(self, ext)

if __name__ == "__main__":
    import sys
    import numpy

    # Publish the library to PyPI.
    if "publish" in sys.argv[-1]:
        os.system("python setup.py sdist upload")
        sys.exit()

    # Choose libraries to link.
    libraries = []
    if os.name == "posix":
        libraries.append("m")

    # Specify the include directories.
    include_dirs = [
        "include",
        numpy.get_include(),
    ]

    # The source files.
    src = [
        "transit/_transit.c",
        "src/quad.cpp",
        "src/driver.cpp",
    ]

    # Set up the extension.
    ext = Extension("transit._transit", sources=src,
                    libraries=libraries, include_dirs=include_dirs)

    # Hackishly inject a constant into builtins to enable importing of the
    # package before the library is built.
    if sys.version_info[0] < 3:
        import __builtin__ as builtins
    else:
        import builtins
    builtins.__TRANSIT_SETUP__ = True
    import transit

    # Execute the setup command.
    desc = open("README.rst").read()
    setup(
        name="transit",
        version=transit.__version__,
        author="Daniel Foreman-Mackey",
        author_email="danfm@nyu.edu",
        packages=["transit"],
        ext_modules=[ext],
        url="http://github.com/dfm/transit",
        license="MIT",
        description="A Python library for computing the light curves of "
                    "transiting planets",
        long_description=desc,
        package_data={"": ["README.rst", "LICENSE", "include/*.h", ]},
        include_package_data=True,
        cmdclass=dict(build_ext=build_ext),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
        ],
    )
