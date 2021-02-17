import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name='pyqchem',
    version='1.0',
    author='Andreas Hauser',
    url='https://github.com/hauser-group/pyQChem',
    license='FreeBSD',
    description='A Python module for scripting with Q-Chem',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    python_requires='>=3.5',
    install_requires=requirements
)
