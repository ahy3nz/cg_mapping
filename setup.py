from setuptools import setup

#####################################
VERSION = "0.0.0"
ISRELEASED = False 
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################


setup(
    name='cg_mapping',
    version=__version__,
    author='Alexander Yang',
    author_email='alexander.h.yang@vanderbilt.edu',
    url='https://github.com/ahy3nz/cg_mapping',
    package_dir={'cg_mapping': 'cg_mapping'},
    license="MIT",
    zip_safe=False,
    keywords='cg_mapping',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
)

