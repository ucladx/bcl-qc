from setuptools import setup

setup( name='bcl-qc',
    version='0.1.0',
    description='Tools to demux and QC data in BCL format',
    url='https://github.com/ucla-dx/bcl-qc',
    author='Cyriac Kandoth',
    author_email='ckandoth@gmail.com',
    license='Apache-2.0',
    packages=['bcl-qc'],
    install_requires=[
        'strsimpy==0.2.1',
    ],
    include_package_data=True,
    zip_safe=False)
