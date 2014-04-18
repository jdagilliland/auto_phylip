from distutils.core import setup

setup(
        name='AutoPHYLIP',
        version='0.1.0',
        author='J. D. A. Gilliland',
        author_email='jdagilliland@gmail.com',
        py_modules=['auto_phylip'],
        scripts=['bin/run_phylip','bin/tab2phy'],
        license='LICENSE.txt',
        description='Utilities to help with using PHYLIP',
        long_description=open('README.txt').read(),
        )
