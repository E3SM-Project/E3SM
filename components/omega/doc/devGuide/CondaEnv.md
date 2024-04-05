(omega-dev-conda-env)=

# Development Conda Environment

As a developer, if you want to build the documentation locally or lint the
code, you will need your own
[conda](https://conda.io/projects/conda/en/latest/index.html) environment with
the necessary packages.

(omega-dev-install-miniforge3)=

## Installing Miniforge3

If you have not already done so, you will install
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3):
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
If you have the space to do so, we recommend installing Miniforge3 somewhere in
our home directory.  Typically, this file system will be faster than on a
scratch space spaces.  Also, a conda environment will be rendered useless if
its contents get purged.  However, conda environments each require several GB
of space, which may quickly use up your disk quota so you will need to weigh
these trade-offs for yourself.

In what follows, we assume you have Miniforge3 installed in
`${HOME}/miniforge3`.

You may wish to skip the step in the Miniforge3 installation where it
wants to modify your `.bashrc`.  If you allow it to, it will add activation of
the `conda` command and the `base` conda environment, and this may not be
something you want.

If you opt not to allow Miniforge3 to update your `.bashrc`, it may be helpful
to add an alias like the following to your `.bashrc`:

```bash
alias conda_base='
source ${HOME}/miniforge3/etc/profile.d/conda.sh
conda activate'
```

Anytime you need `conda` command, you can first run `conda_base`to get it.


(omega-dev-create-conda-env)=

## Creating the Conda Environment

With Miniforge3 installed and the `conda` command available, you can create the
Omega development environment, starting from the root
of an Omega branch, with:

```bash
cd components/omega/
conda create -n omega_dev --file dev-conda.txt
conda activate omega_dev
```

(omega-dev-update-conda-env)=

## Updating the Conda Environment

If you need to update the `omega_dev` conda environment, the easiest way to
do that is just to recreate it by following the same instructions as above.
You will be prompted about whether you want to delete the old environment
unless you include the `-y` flag (to say "yes" to all prompts).
