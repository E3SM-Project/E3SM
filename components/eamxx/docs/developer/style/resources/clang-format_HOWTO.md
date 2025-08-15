# How-to Guide: `clang-format`

This guide is for developers who wish to apply `clang-format` on their chosen
development machine, whether that be their personal machine or a multi-user
cluster.
In this guide, we will describe how to configure/run `clang-format` on EAMxx
code and also how to install it if necessary.

## Configure and Run `clang-format`

Running `clang-format` according to the defined EAMxx standard ***can*** be
done using only command line arguments; however, the command is quite long.
The easier route is to reference the configuration file
(`$EAMXX_ROOT/.clang-format`).
In this case the command is

```bash {.copy}
clang-format [-i] --style="file:${EAMXX_ROOT}/.clang-format" <path/to/c++/source/file(s)>
```

where the `-i` (`--in-place`) argument controls whether the formatting edits
are conducted and change the source file(s) ("in place") or the required edits
are printed to `stdout` (flag omitted).
Also, note that the `--style` flag can also fully define the options without
the configuration file as follows, but this is not recommended as the
configuration could change before this guide is updated to reflect as such.

<!-- markdownlint-disable MD046 MD013 -->
??? Danger "Terminal One-(long-)liner"

    ```bash {.copy}
    clang-format [-i] --style="{BasedOnStyle: llvm, ColumnLimit: 100, AlignConsecutiveAssignments: true, AlignConsecutiveBitFields: true, AlignConsecutiveMacros: true, AlignEscapedNewlines: true, AlignTrailingComments: true}" <path/to/c++/source/file(s)>
    ```
<!-- markdownlint-enable MD046 MD013 -->

## Installing `clang-format`

On a personal machine, which we will consider to be one for which you have
`sudo` privileges, installation can be conducted via package manager or by
building `llvm v14` from scratch.
If you are a non-admin user of a multi-user cluster, HPC platform, etc., then
it is likely to be an ***easy*** process, though potentially not immediate.

<!-- markdownlint-disable MD046 -->
??? Note "If you require or already have multiple versions of `clang-format` installed"

    Note that, depending on your requirements, this could be placed in your shell
    config file (`.bashrc`, `.zshrc`), or if only needed for a shell session, it
    can be done directly from the command line.

    The other option is to add a versioned symbolic link to your `PATH`.
    This is sometimes included in the package's `bin/` directory by default and,
    if not, can be added there after of placed somewhere that is already on your
    `PATH`.

    ```c++
    $ cd /opt/homebrew/opt/llvm@14/bin
    $ ls clang-format*
    clang-format
    // no versioned binary so we create symlink
    $ ln -s ./clang-format ./clang-format-14
    $ which clang-format
    /opt/homebrew/opt/llvm@14/bin/clang-format-14
    ```

    **OR**

    ```c++
    $ cd /opt/homebrew/opt/llvm@14/bin
    $ ls clang-format*
    clang-format
    // no versioned binary so we create symlink in a directory already on PATH
    $ echo $PATH
    /usr/bin:/usr/local/bin
    $ ln -s ./clang-format /usr/local/bin/clang-format-14
    $ which clang-format
    /usr/local/bin/clang-format-14
    ```
<!-- markdownlint-enable MD046 -->

=== "Personal Machine"

    <!-- markdownlint-disable MD046 -->
    === "Mac Users (Homebrew Package Manager)"

        For Mac users, the [Homebrew](https://brew.sh)
        package manager (`brew`) is the quickest and most straightforward way
        to get `clang-format` installed.
        Since `clang-format` v14 is not available to install directly from
        Homebrew, we install the entire LLVM package at version 14, and this
        is as simple as

        ```bash {.copy}
        brew install llvm@14
        ```

        It it likely that Homebrew will issue a message about not linking the
        `clang`-related tools by default, so next we add the binary to our `PATH`.

        ```bash {.copy}
        $ export PATH="/opt/homebrew/opt/llvm@14/bin/clang-format:$PATH"
        # Note: this is the default location for a recent Mac running Apple silicon.
        # It may be different on another system.
        # You can confirm where yours is installed via 'brew info llvm@14'
        $ which clang-format
        /opt/homebrew/opt/llvm@14/bin/clang-format
        ```

        Note also, that if your system has multiple version of `clang-format` installed,
        it may be preferable to instead set a versioned alias to `clang-format`
        (`clang-format-14`) as in

        ```c++
        // create a shell-alias
        alias clang-format-14="/opt/homebrew/opt/llvm@14/bin/clang-format"
        ```

    === "Linux Users (Package Manager)"

        Given the many flavors of Linux, it is difficult to generalize, but there
        is a high probability the proper version of `clang-format` or `llvm` is
        provided by the built-in package manager.
        The commands will differ based on your Linux distribution, but using the
        Debian/Ubuntu `apt` syntax, it could be accomplished via something like

        ```bash {.copy}
        $ apt search clang-format
        [...]
        clang-format-14/...
        $ apt install clang-format-14
        ```

    === "Build from Source"

        If you do not succeed in the above, `clang-format` can also be fully built
        from the [LLVM Compiler Infrastructure](https://github.com/llvm/llvm-project).
        It will begin with something like

        ```bash {.copy}
        git clone git@github.com:llvm/llvm-project.git
        git checkout llvmorg-14.0.6 # version tag
        ```

        Also, if you only need `clang-format` and not any of the other tooling,
        it will build faster/smaller if you use the CMake flag
        `-DLLVM_ENABLE_PROJECTS="clang"` to only build `clang` and it's friends,
        rather than all of `llvm`.
        And finally, the README for [LLVM version 14.0.6](https://github.com/llvm/llvm-project/tree/llvmorg-14.0.6)
        is far more comprehensive that the one for the latest version, and it contains
        instructions specific to that build.

=== "Multi-user System"

    In many cases `llvm`, `clang`, or `clang-format` will be available as a module,
    though whether version 14 is an option could be less likely.
    In the optimistic case, it could be as simple as (using Lmod syntax)
    
    ```bash {.copy}
    $ module avail llvm # [clang, clang-format]
    [... list ]
    $ module load llvm/14.0.6
    ```

    If it is not available, you will probably need to reach out to a system
    administrator to get an official version installed.[^but-conda]
<!-- markdownlint-enable MD046 -->

---

<!-- markdownlint-disable MD046 -->
??? Tip "Unnecessary but Convenient Workflow Customization (`direnv`)"

    If you'd like to add a layer of automation/complexity to ensure you only use
    `clang-format v14` on `EAMxx` and/or want to use a newer version on the rest
    of your system, there is a very handy terminal tool called
    [direnv](https://direnv.net/)
    (`brew install direnv`) that allows you to automatically load and unload
    environment variables based on `$PWD` using a `.envrc` file.
    As an example, here's my `.envrc` that adds `clang-format v14` to the path
    when I'm working on `EAMxx`.

    ```bash {.copy}
    PATH_add /opt/homebrew/opt/llvm@14/bin/clang-format

    # also, since I often forget to load a python environment that's required for
    # running ctest, this creates or loads a python 3 virtual environment with numpy
    layout python3
    pip install --upgrade pip
    # the upgrade isn't strictly necessary but trades a little extra setup on
    # the front end to avoid pip endlessly reminding you to update
    pip install numpy
    ```

    This file can be placed in the top-level `EAMxx` directory, and running
    `direnv allow` enables the functionality.
    Namely, executing `cd EAMxx` loads the defined environment that stays loaded
    in any subdirectories, and resets the standard environment when exiting to a
    directory above/outside of `EAMxx`.

    For the `conda` fans, this tool can also be used to auto-load a
    pre-configured `conda` environment since the `.envrc` is essentially a bash
    script with a few bells and whistles tacked on.
<!-- markdownlint-enable MD046 -->

<!-- markdownlint-disable MD053 -->
[^but-conda]: There are rumors of using `conda` creatively to do a user-install,
but that is not an option we support or suggest.
<!-- markdownlint-enable MD053 -->
