# EAMxx Code Formatting Standards

To enforce consistent code format throughout EAMxx, we make use of an
autoformatting workflow, carried out via Github Actions in the E3SM repository.

- The tool we employ is
  [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html),
  and the version we have chosen is v14.[^v14]
<!-- FIXME: Nail this down when figured out -->
- The standard we maintain is largely identical to the
  [LLVM Coding Standards](https://llvm.org/docs/CodingStandards.html),
  and a list of the handful of customizations of this convention are
  enumerated in the configuration file [`$EAMXX_ROOT/.clang-format`](https://github.com/E3SM-Project/E3SM/blob/master/components/eamxx/.clang-format).
      - See this [How-to Guide](resources/clang-format_HOWTO.md)
        for additional details on how to configure `clang-format` on your chosen
        development machine.

## Automated Workflow Summary

- The `eamxx-format` workflow runs automatically and passes or fails based on
  adherence to our formatting standard.
- The workflow is triggered by any Pull Request (PR) that modifies EAMxx code.
      - It is also triggered by other PR-related triggers, such as pushing
        changes, or converting from **draft** to **ready**.
- All code modified by a PR must be formatted prior to **merging**.
- It is not necessary for your code to be formatted upon ***opening*** a
  Pull Request, but feel free to `clang-format` as you develop if that is
  your preference.
      - The one situation for which opening a pre-formatted PR may not be
        preferred is if the file has never previously been `clang-format`-ed
        and requires a large number of changes.
            - I.e., touching the majority of lines in a file for format-only
              changes will make it difficult for a reviewer to determine which
              lines were changed for substantive reasons.
            - In this case, please refrain from formatting the code prior to
              opening the PR, and, instead, run `clang-format` once the PR is
              approved to be merged and make that the final commit.
- As of now, the `eamxx-format` workflow only considers files that are edited
  by the Pull Request.
- In addition to the pass/fail status, the workflow, provides the path to any
  files that caused the failure.
      - This information is found on the ***Summary*** page for any failed
        `eamxx-format` ***Job***.[^huh-where]

[^v14]: It turns out that this is important because there really are
differences in behavior across versions.
[^huh-where]: To get to this summary, select the ***Checks*** tab at the top of
the PR page and select the `eamxx-format` workflow from the left sidebar.
The summary is in the main pane of the page with the title
**clang-format-linter summary**.[^also]
[^also]: Note that this can also be accessed the long way around by following the
***Actions*** link at the top of the E3SM repository page;
select `eamxx-format` from the ***All workflows*** section of the ***Actions***
sidebar; then choose the most recent run that is associated with your PR,
which should be near the top of the list.
