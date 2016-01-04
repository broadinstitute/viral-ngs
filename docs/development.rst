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

    $ pyenv install 2.7.10 3.4.3 3.5.0

3. Change directory to the `viral-ngs` git checkout.
4. Set the local pyenv versions:

::

    $ pyenv local 3.5.0 3.4.3 2.7.10

5. Install tox and tox-pyenv

::

    $ pip3 install tox tox-pyenv

6. Run tox to run tests and check building docs.

::

    $ tox

Testing easy_deploy with Vagrant and Ansible
--------------------------------------------

By default, easy_deploy installs a version of viral-ngs by cloning the github
repository, which makes it difficult to test local changes to the playbook.
Testing the easy_deploy setup locally can be done using Vagrant and the Ansible
playbook with some custom commands. To do so, we need to take the following
steps:

1. `Install Vagrant <https://docs.vagrantup.com/v2/installation/>`_
2. `Install Ansible <http://docs.ansible.com/ansible/intro_installation.html>`_
3. Change to the easy_deploy directory

::

   $ cd easy_deploy

4. Run the easy deploy script to bootstrap the VM and answer the prompts

::

   $ ./run.sh

5. From now on, run Ansible playbook manually for provisioning:
   http://docs.ansible.com/ansible/guide_vagrant.html#running-ansible-manually.

6. Add `deploy=sync` or `deploy=archive` to the `--extra-vars` of the Ansible
   playbook command

::

   $ ansible-playbook ... --extra-vars=deploy=sync

To test playbook changes in a local repository, the choice of `deploy=sync` or
`deploy=archive` changes the strategy used to "deploy" viral-ngs into the
Vagrant VM.

The `deploy=sync` option creates a symlink to the synced folder containing the
root of your viral-ngs git repo. Therefore, any tool installation or other
changes will be reflected on both the host and the guest machine. This is
desirable for fast iteration of changes, but makes it difficult to isolate the
host's viral-ngs installation from the installation on the guest VM.

The `deploy=archive` option performs a `git archive` on the host's viral-ngs
repository and untars it into the project directory. This is a clean install
each time, which can be time consuming due to the need to reinstall all
dependencies, and only the current HEAD commit will be reflected in the guest.
No uncommitted/dirty changes will be picked up using this method. This is more
ideally suited for a finalized clean test of the playbook.
