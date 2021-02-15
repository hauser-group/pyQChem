import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
    
install_requires = [
    'numpy>=1.10.4',
    'matplotlib>=2.0.0',
]

setuptools.setup(
    name='pyqchem',
    version='0.0.2',
    author='Andreas Hauser',
    url='https://github.com/hauser-group/pyQChem',
    license='FreeBSD',
    description='A Python module for scripting with Q-Chem',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    python_requires='>=3.6'
)
