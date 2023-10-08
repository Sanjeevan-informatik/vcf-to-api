from setuptools import setup

setup(
    name='vcftoapi',
    packages=['vcftoapi'],
    include_package_data=True,
    install_requires=[
        'flask',
    ],
)