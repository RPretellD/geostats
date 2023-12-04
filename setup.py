from setuptools import setup, find_packages, Extension

with open("README.md","r", encoding = 'utf-8') as fp:
	readme = fp.read()

setup(
	name="geostats",
	version="0.1.0",
	description="A suite of geostatistical tools.",
	author="A. Renmin Pretell Ductram",
	author_email='rpretell@unr.edu',
	url="https://github.com/RPretellD/geostats",
    long_description=readme,
    
    packages=find_packages(),
	include_package_data=True,
	
    install_requires=["numpy","Cython"],

	license='MIT',
	keywords='geostats',
	classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
	]
)