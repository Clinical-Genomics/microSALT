import codecs

from setuptools import setup

def parse_reqs(req_path='./requirements.txt'):
    """Recursively parse requirements from nested pip files."""
    install_requires = []
    with codecs.open(req_path, 'r') as handle:
        # remove comments and empty lines
        lines = (line.strip() for line in handle
                 if line.strip() and not line.startswith('#'))
        for line in lines:
            # check for nested requirements files
            if line.startswith('-r'):
                # recursively call this function
                install_requires += parse_reqs(req_path=line[3:])
            else:
                # add the line as a new requirement
                install_requires.append(line)
    return install_requires

setup(
    name="locioser",
    version="1.0.0",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=parse_reqs(),  

    entry_points={
        'console_scripts': [
            'locioser=core:main'  
        ]
    },
)
