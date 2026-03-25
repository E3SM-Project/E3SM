# Composite action to update merge commit when in a pull request event

When re-running a workflow that was run as part of PR testing (i.e., the
github event was `pull_request`), github uses _the same_ merge commit,
in an attempt to re-run _the same_ job as the previous run.
However, this is often NOT what the user wants. The user wants to re-run
the workflow on the current merge commit, to include any modification
(including bug fixes) that were merged into master in the meantime.

This action allows to check the state of the merge commit compared to
the target branch (usually master), to ensure that the merge commit
is "ahead" of the target branch (meaning that last commit in the
target branch is an ancestor of the merge commit). If not, proceed
to merge the target branch into HEAD.
