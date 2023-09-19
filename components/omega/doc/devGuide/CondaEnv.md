(omega-dev-conda-env)=

# Development Conda Environment

As a developer, if you want to build the documentation locally or lint the
code, you will need your own
[conda](https://conda.io/projects/conda/en/latest/index.html) environment with
the necessary packages.

(omega-dev-install-mambaforge)=

## Installing Mambaforge

If you have not already done so, you will install
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge), preferably
somewhere in your home directory but possibly somewhere in your scratch space.
In what follows, we assume you have Mambaforge installed in
`${HOME}/mambaforge`.

You may wish to skip the step in the Mambaforge installation where it
wants to modify your `.bashrc`.  If you allow it to, it will add
activation of the `mamba` and `conda` commands and the `base` conda
environment, and this may not be something you want.

If you opt not to allow Mambforge to update your `.bashrc`, it may be helpful
to add an alias like the following to your `.bashrc`:

```bash
alias mamba_base='
source ${HOME}/mambaforge/etc/profile.d/conda.sh
source ${HOME}/mambaforge/etc/profile.d/mamba.sh
mamba activate'
```

Anytime you need `conda` and `mamba` commands, you can first run `mamba_base`
to get them.


(omega-dev-create-conda-env)=

## Creating the Conda Environment

With Mambaforge installed and the `conda` and `mamba` commands activated,
you can create the OMEGA development environment, starting from the root
of an OMEGA branch, with:

```bash
cd components/omega/
mamba create -n omega_dev --file dev-conda.txt
conda activate omega_dev
```
