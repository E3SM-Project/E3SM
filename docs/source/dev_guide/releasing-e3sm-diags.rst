.. _prepare-release:

How to Prepare a Release
========================

In this guide, we'll cover:

* Bumping the Version
* Releasing on Github
* Releasing on Anaconda
* Creating a New Version of the Documentation

Bumping the Version
-------------------

1. Checkout the latest ``main``.
2. Checkout a branch with the name of the version.

    ::

        # Prepend "v" to <version>
        # For release candidates, append "rc" to <version>
        git checkout -b bump/v<version>
        git push --set-upstream origin bump/v<version>

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
        - conda-env/prod.yml:22 - e3sm_diags=2.5.0
        + conda-env/prod.yml:22 - e3sm_diags=2.6.0
        - tbump.toml:5 current = "2.5.0"
        + tbump.toml:5 current = "2.6.0"
        => Would run these git commands
        $ git add --update
        $ git commit --message Bump to 2.6.0
        $ git push origin bump/v2.6.0
        :: Looking good? (y/N)
        >

4. Create a pull request to the main repo and merge it.

.. _github-release:

Releasing on GitHub
-------------------

1. Draft a new release on the `releases page <https://github.com/E3SM-Project/e3sm_diags/releases>`_.
2. Set `Tag version` to ``v<version>``, **including the "v"**. `@Target` should be ``main``.
3. Set `Release title` to ``v<version>``, **including the "v"**.
4. Use `Describe this release` to summarize the changelog.

   * You can scroll through `e3sm-diags commits <https://github.com/E3SM-Project/e3sm_diags/commits/main>`_ for a list of changes.

5. If this version is a release candidate (``<version>`` appended with ``rc``), checkmark `This is a pre-release`.
6. Click `Publish release`.
7. CI/CD release workflow is automatically triggered.

Releasing on Anaconda
---------------------

1. Be sure to have already completed :ref:`Releasing on GitHub <github-release>`.
2. Follow the steps `on the feedstock page <https://github.com/conda-forge/e3sm_diags-feedstock#updating-e3sm_diags-feedstock>`_ for updating the ``e3sm_diags-feedstock``.
3. The package will be released on ``conda-forge``.

Creating a New Version of the Documentation
-------------------------------------------

1. Be sure to have already completed :ref:`Creating A Release On GitHub <github-release>`. This triggers the CI/CD workflow that handles publishing documentation versions.
2. Wait until the CI/CD build is successful. You can view all workflows at `All Workflows <https://github.com/E3SM-Project/e3sm_diags/actions>`_.
3. Changes will be available on the `e3sm_diags documentation page <https://e3sm-project.github.io/e3sm_diags/>`_.
