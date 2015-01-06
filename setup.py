from setuptools import setup, find_packages
import util.version

def read_file(fname):
    with open(fname, 'rt') as inf:
        return list(x.rstrip('\n\r') for x in inf)

setup(
    name='viral_ngs',
    version=util.version.get_version(),
    url='http://github.com/broadinstitute/viral-ngs/',
    license='BSD-style Software License',
    author='Broad Viral Genomics / Sabeti Lab',
    tests_require=['nose'],
    install_requires=read_file('requirements.txt'),
    author_email='dpark@broadinstitute.org',
    description='NGS analysis pipelines for viral projects',
    #long_description='''x''',
    packages=['viral_ngs'],
    include_package_data=True,
    platforms='any',
    test_suite='nose.collector',
    extras_require={
        'pipeline': ['snakemake'],
    }
)
