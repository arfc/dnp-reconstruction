from setuptools import setup, find_packages

setup(
    name="mosden",
    version="0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'mosden = mosden.__main__:main'
        ]
    },
)
