Developing viral-ngs
====================

Testing with tox and pyenv
--------------------------

Using pyenv with tox can simplify testing locally against multiple different Python
versions and dependency environments. To setup your development environment for
testing, you first need to perform the following steps:

1. Install pyenv for your OS https://github.com/yyuu/pyenv#installation
2. Install the appropriate Python versions:

::

    pyenv install 2.7.10 3.4.3 3.5.0

3. Change directory to the `viral-ngs` git checkout.
4. Set the local pyenv versions:

::

    pyenv local 3.5.0 3.4.3 2.7.10

5. Install tox and tox-pyenv

::

    pip3 install tox tox-pyenv

6. Run tox to run tests and check building docs.

::

    tox

