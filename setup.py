from setuptools import setup, find_packages


with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="crystalmanipulator",
    version="0.1",
    packages=find_packages(),
    install_requires=required,
    include_package_data=True,
    package_data={
        '': ['BASIC_DATA'],
    },
    entry_points={
        'console_scripts': [
            'crystalmanipulator=crystalmanipulator.main:main'
        ]
    },    
    author="Schrodinger's Cat",
    author_email="gmy721212@163.com",
    description="CrystalManipulator",
    keywords="Crystal Manipulator",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

