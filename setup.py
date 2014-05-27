from distutils.core import setup

setup(
        name='AutoPHYLIP',
        version='0.1.1',
        author='J. D. A. Gilliland',
        author_email='jdagilliland@gmail.com',
        py_modules=['auto_phylip'],
        scripts=[
            'bin/run_phylip',
            'bin/tab2phy',
            'bin/run_seqboot',
            'bin/run_consense',
            'bin/cleanup_consense',
            ],
        license='LICENSE.txt',
        description='Utilities to help with using PHYLIP',
        long_description=open('README.txt').read(),
        )
