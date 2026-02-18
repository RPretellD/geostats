from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.errors import CCompilerError, DistutilsExecError, DistutilsPlatformError
from Cython.Build import cythonize
import numpy

with open("README.md","r", encoding = 'utf-8') as fp:
	readme = fp.read()

# Define the Cython extensions
extensions = [

    Extension(
        name = "geostats.build_Mrho_cy",
        sources = [ "geostats/build_Mrho_cy.pyx",
        ],
        include_dirs = [numpy.get_include()],
    ),

    Extension(
        name = "geostats.geostats_tools_cy",
        sources = [ "geostats/geostats_tools_cy.pyx",
        ],
        include_dirs = [numpy.get_include()],
    ),
]

class OptionalBuildExt(build_ext):

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            self.extensions = []

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError):
            pass

# Call setup
setup(
    name="geostats",
    version="1.0.0",
    license='MIT',
    
    description="A suite of geostatistical tools.",
    author="A. Renmin Pretell Ductram",
    author_email='rpretell@unr.edu',
    url="https://github.com/RPretellD/geostats",
    
    long_description_content_type="text/markdown",
    long_description=readme,
    
    packages=find_packages(),

	include_package_data=True,
    ext_modules=cythonize(extensions),
    cmdclass={"build_ext": OptionalBuildExt},
    python_requires	= ">=3.7",
    install_requires=[
            "numpy",
            "Cython",
            "vincenty",
    ],

	keywords='geostats',
	classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
	],    
)
