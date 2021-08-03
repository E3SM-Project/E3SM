.. _prepare-release:

How to Prepare a Release
========================

In this guide, we'll cover:

* Bumping the Version
* Releasing on Github
* Releasing on Anaconda
* Creating a New Version of the Documentation
* Building and Releasing the Docker Image

Bumping the Version
-------------------

1. Checkout the latest ``master``.
2. Checkout a branch with the name of the version.

    ::

        # Prepend "v" to <version>
        # For release candidates, append "rc" to <version>
        git checkout -b v<version>
        git push --set-upstream origin v<version>

3. Bump version using tbump.

    ::

        # Exclude "v" and <version> should match step 2
        # --no-tag is required since tagging is handled in "Releasing on GitHub"
        tbump <version> --no-tag

        :: Bumping from 2.5.0 to 2.6.0
        => Would patch these files
        - setup.py:113 version="2.5.0",
        + setup.py:113 version="2.6.0",
        - acme_diags/__init__.py:4 __version__ = "v2.5.0"
        + acme_diags/__init__.py:4 __version__ = "v2.6.0"
        - conda/meta.yaml:2 {% set version = "2.5.0" %}
        + conda/meta.yaml:2 {% set version = "2.6.0" %}
        - conda/e3sm_diags_env.yml:22 - e3sm_diags=2.5.0
        + conda/e3sm_diags_env.yml:22 - e3sm_diags=2.6.0
        - Dockerfile:8 LABEL version="2.5.0"
        + Dockerfile:8 LABEL version="2.6.0"
        - tbump.toml:5 current = "2.5.0"
        + tbump.toml:5 current = "2.6.0"
        => Would run these git commands
        $ git add --update
        $ git commit --message Bump to 2.6.0
        $ git push origin v2.6.0
        :: Looking good? (y/N)
        >

4. Create a pull request to the main repo and merge it.

.. _github-release:

Releasing on GitHub
-------------------

1. Draft a new release on the `releases page <https://github.com/E3SM-Project/e3sm_diags/releases>`_.
2. Set `Tag version` to ``v<version>``, **including the "v"**. `@Target` should be ``master``.
3. Set `Release title` to ``v<version>``, **including the "v"**.
4. Use `Describe this release` to summarize the changelog.

   * You can scroll through `e3sm-diags commits <https://github.com/E3SM-Project/e3sm_diags/commits/master>`_ for a list of changes.

5. If this version is a release candidate (``<version>`` appended with ``rc``), checkmark `This is a pre-release`.
6. Click `Publish release`.
7. CI/CD release workflow is automatically triggered.

Releasing On Anaconda
---------------------

1. Be sure to have already completed :ref:`Releasing on GitHub <github-release>`.
2. Follow the steps `on the feedstock page <https://github.com/conda-forge/e3sm_diags-feedstock#updating-e3sm_diags-feedstock>`_ for updating the ``e3sm_diags-feedstock``.
3. The package will be released on ``conda-forge``.

Creating a New Version of the Documentation
-------------------------------------------

1. Be sure to have already completed :ref:`Creating A Release On GitHub <github-release>`. This triggers the CI/CD workflow that handles publishing documentation versions.
2. Wait until the CI/CD build is successful. You can view all workflows at `All Workflows <https://github.com/E3SM-Project/e3sm_diags/actions>`_.
3. Changes will be available on the `e3sm_diags documentation page <https://e3sm-project.github.io/e3sm_diags/>`_.

How To Build and Release The Docker Image
-----------------------------------------

A Docker image of ``e3sm_diags`` needs to be created and released as well.
This Docker image can be ran as a container via Docker, Shifter, or Singularity.

We'll build the image, test it, and then release it.

Prerequisites
^^^^^^^^^^^^^

1. Please make a Docker ID if you haven't done so already.
This is needed to release and upload the image.

2. Also make sure that you have access to the `e3sm Dockerhub <https://hub.docker.com/u/e3sm>`_ ,
and specifically the e3sm_diags repo there. If you don't, you'll see an error when you run
``docker push`` later on in this guide.
Email Jill Zhang (zhang40@llnl.gov) or Rob Jacob (jacob@anl.gov) for access.


Building
^^^^^^^^

3. Set an environmental variable, ``E3SM_DIAGS_VERSION``, to the version that you're releasing.

    ::

        export E3SM_DIAGS_VERSION=v1.5.0

A Temporary Diversion
"""""""""""""""""""""

4. When installing the software, a user needs to do ``pip install --user .``
instead of the traditional ``python setup.py install``.
It's the way Anaconda recommends creating packages.
This is *currently* causing issues when building the Docker image.
Due to this, open ``setup.py`` and change the ``INSTALL_PATH`` to be ``os.path.join(sys.prefix, 'share/e3sm_diags/')``.

    .. code-block:: python

        # INSTALL_PATH = 'share/e3sm_diags/'
        INSTALL_PATH = os.path.join(sys.prefix, 'share/e3sm_diags/')


5. Open the ``Dockerfile`` and change any instance of ``pip install --user .`` to ``python setup.py install``.

    ::

        RUN conda env update -n base --file conda/e3sm_diags_env_dev.yml && \
                conda clean --all -y && \
                source activate base && \
                # pip install --user . && \
                python setup.py install && \
                rm -r build/

Back to Building the Image
""""""""""""""""""""""""""

6. Go to the root of the project, where the ``Dockerfile`` is located and run the command below.
This builds the image and adds two tags, one titled ``latest`` and one based on the version you're releasing.
By prefixing the tag with ``e3sm/``, it'll upload it to the
`e3sm Dockerhub <https://hub.docker.com/u/e3sm>`_,
which we'll do in forthcoming steps.

When Docker builds an image, it sends all of the data in the current working directory as the build context.
So if the current directory has a lot of data (like sample runs, large nc files, etc),
remove them before continuing.
Check the size of the current directory with ``du -sh .``.

    ::

        docker build . -t e3sm/e3sm_diags:latest -t e3sm/e3sm_diags:$E3SM_DIAGS_VERSION


7. View the Docker images you have. You should see the images you've made, based on the tags.

    ::

        docker images

You should see something like this:

    ::

        REPOSITORY               TAG                 IMAGE ID            CREATED             SIZE
        e3sm/e3sm_diags          latest              bc7f93375025        6 minutes ago       3.57GB
        e3sm/e3sm_diags          v1.5.0              bc7f93375025        6 minutes ago       3.57GB
        continuumio/miniconda    4.5.4               16e4fbac86ce        7 weeks ago         544MB
        hello-world              latest              e38bc07ac18e        5 months ago        1.85kB

Testing
"""""""

8. Go to the folder with the system tests.

    ::

        cd tests/system/


9. ``wget`` or ``curl`` the script to run the image.
When you actually run an image, it's called a **container**.

    ::

        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py

        # Or use this:
        curl -O https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py


10. Run the tests. Check the terminal and
results after each run to ensure that everything was created without errors.

    ::

        python e3sm_diags_container.py --docker -p all_sets.py -d all_sets.cfg


11. If you do find an error, it could be with the script ``e3sm_diags_container.py`` or with ``e3sm_diags`` itself.
Please fix this. You might need to delete the release, or release a bug-fix version.

Releasing
"""""""""

12. Push both of the images, one with the ``latest`` tag and the other with the version you're releasing.

::

    docker push e3sm/e3sm_diags:latest
    docker push e3sm/e3sm_diags:$E3SM_DIAGS_VERSION


13. Congratulations, you're done! You can go home/nap for the day, I won't tell.

Optional: Cleanup
"""""""""""""""""

* These images can take up a fair amount of space on your machine, since each is around 4GB.
  Here are some ways to manage them.

  * View all of the images you have with ``docker images``.
    You can remove an image by the image id.
    The ``--force`` option is also supported.

    ::

        docker rmi <image_id>

  * Run the command below once in a while to remove unused data.
    This includes any intermediate or broken images/container.

    ::

        docker system prune
