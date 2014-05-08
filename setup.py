#!/usr/bin/env python
"""Implementations of algorithms related to the work of Miroslav Fiedler.

"""

DOCLINES = __doc__.split('\n')

# This setup script is written according to
# http://docs.python.org/2/distutils/setupscript.html
#
# It is meant to be installed through github using pip.

from distutils.core import setup

setup(
        name='fiedlerology',
        version='0.1',
        description=DOCLINES[0],
        author='alex',
        url='https://github.com/argriffing/fiedlerology/',
        download_url='https://github.com/argriffing/fiedlerology/',
        packages=['fiedlerology'],
        test_suite='nose.collector',
        package_data={'fiedlerology' : ['tests/test_*.py']},
        )


