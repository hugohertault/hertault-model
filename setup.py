```python
from setuptools import setup, find_packages

setup(
    name="hertault-model",
    version="1.0.0",
    author="Hugo Hertault",
    description="Unified dark sector model with environmental phase transition",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/hugohertault/hertault-model",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0", 
        "matplotlib>=3.3.0",
        "pandas>=1.3.0"
    ],
)
```
