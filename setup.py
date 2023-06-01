from setuptools import setup, find_packages

setup(
        name='madx',
        version='0.0.0',
        description='Methodical Accelerator Design MAD-X',
        author='CERN',
        author_email='riccardo.de.maria@cern.ch',
        url='https://github.com/MethodicalAcceleratorDesign/mad-x',
        packages=find_packages(),
        install_requires=['numpy','matplotlib','cpymad'],
)
