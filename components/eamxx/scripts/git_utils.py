"""
Git-related Utilities
"""

from utils import expect, run_cmd, run_cmd_no_fail

import os

###############################################################################
def get_current_branch(repo=None):
###############################################################################
    """
    Return the name of the current branch for a repository
    If in detached HEAD state, returns None
    """

    stat, output, err = run_cmd("git rev-parse --abbrev-ref HEAD", from_dir=repo)
    expect (stat==0, "Error! The command 'git rev-parse --abbrev-ref HEAD' failed with error: {}".format(err))

    return None if output=="HEAD" else output

###############################################################################
def get_current_commit(short=False, repo=None, tag=False, commit="HEAD"):
###############################################################################
    """
    Return the sha1 of the current HEAD commit

    >>> get_current_commit() is not None
    True
    """
    # This is for testing only
    if "SCREAM_FAKE_GIT_HEAD" in os.environ:
        return os.environ["SCREAM_FAKE_GIT_HEAD"]

    if tag:
        rc, output, err = run_cmd("git describe --tags $(git log -n1 --pretty='%h')", from_dir=repo)
    else:
        rc, output, err = run_cmd("git rev-parse {} {}".format("--short" if short else "", commit), from_dir=repo)

    if rc != 0:
        print("Warning: getting current commit {} failed with error: {}".format(commit, err))

    return output if rc == 0 else None

###############################################################################
def get_current_head(repo=None):
###############################################################################
    """
    Return current head, preferring branch name if possible
    """
    branch = get_current_branch(repo=repo)
    if not branch:
        return get_current_commit(repo=repo)
    else:
        return branch

###############################################################################
def git_refs_difference (cmp_ref, head="HEAD", repo=None):
###############################################################################
    """
    Return the difference in commits between cmp_ref and head.
    In particular, it returns two numbers: the number of commits
    in cmp_ref that are not in head, and the number of commits in head
    that are not in cmp_ref. The former is how much head is behind cmp_ref,
    while the latter is how much head is ahead of cmp_ref.
    """
    if "SCREAM_FAKE_GIT_HEAD" in os.environ:
        expect("SCREAM_FAKE_AHEAD" in os.environ,
               "git_refs_difference cannot be used with SCREAM_FAKE_GIT_HEAD and without SCREAM_FAKE_AHEAD")
        return 0, 0 if cmp_ref == head else int(os.environ["SCREAM_FAKE_AHEAD"])

    cmd = "git rev-list --left-right --count {}...{}".format(cmp_ref, head)
    out = run_cmd_no_fail("{}".format(cmd), from_dir=repo)

    behind_ahead = out.split()
    expect(len(behind_ahead)==2, "Error! Something went wrong when running {}".format(cmd))
    behind, ahead = int(behind_ahead[0]), int(behind_ahead[1])

    return behind, ahead

###############################################################################
def is_repo_clean(repo=None, silent=False):
###############################################################################
    rc, output, _ = run_cmd("git status --porcelain --untracked-files=no", combine_output=True, from_dir=repo)
    if (rc != 0 or output != "") and not silent:
        print("Warning: repo is not clean: {}".format(output))

    return rc == 0 and output == ""

###############################################################################
def get_common_ancestor(other_head, head="HEAD", repo=None):
###############################################################################
    """
    Returns None on error.
    """
    rc, output, _ = run_cmd("git merge-base {} {}".format(other_head, head), from_dir=repo)
    return output if rc == 0 else None

###############################################################################
def update_submodules(repo=None):
###############################################################################
    """
    Updates submodules
    """
    rc, _, errput = run_cmd("git submodule update --init --recursive", from_dir=repo)
    if rc != 0:
        print("Warning: normal submodule update failed:\n{}".format(errput))

        # Trying again without fetching because fetch can fail without ssh keys set up
        run_cmd_no_fail("git submodule update --no-fetch --init --recursive", from_dir=repo)

###############################################################################
def merge_git_ref(git_ref, repo=None, verbose=False, dry_run=False):
###############################################################################
    """
    Merge given git ref into the current branch, and updates submodules
    """

    # Even thoguh it can allow some extra corner cases (dirty repo, but ahead of git_ref),
    # this check is mostly for debugging purposes, as it will inform that no merge occurred
    out = get_common_ancestor(git_ref)
    if out == get_current_commit(commit=git_ref):
        if verbose:
            print ("Merge of '{}' not necessary. Current HEAD is already ahead.".format(git_ref))
        return

    merge_cmd = "git merge {0} -m 'Automatic merge of {0}'".format(git_ref)
    if dry_run:
        print("Would run: {}".format(merge_cmd))
    else:
        expect(is_repo_clean(repo=repo), "Cannot merge ref '{}'. The repo is not clean.".format(git_ref))
        run_cmd_no_fail(merge_cmd, from_dir=repo)
        update_submodules(repo=repo)
        expect(is_repo_clean(repo=repo), "Something went wrong while performing the merge of '{}'".format(git_ref))
        if verbose:
            print ("git ref {} successfully merged.".format(git_ref))
            print_last_commit()

###############################################################################
def print_last_commit(git_ref=None, repo=None, dry_run=False):
###############################################################################
    """
    Prints a one-liner of the last commit
    """
    if dry_run:
        print("Last commit on ref '{}'".format(git_ref))
    elif "SCREAM_FAKE_GIT_HEAD" in os.environ:
        print("Last commit on ref '{}'".format(os.environ["SCREAM_FAKE_GIT_HEAD"]))
    else:
        git_ref = get_current_head(repo) if git_ref is None else git_ref
        last_commit = run_cmd_no_fail("git log --oneline -1 {}".format(git_ref), from_dir=repo)
        print("Last commit on ref '{}': {}".format(git_ref, last_commit))

###############################################################################
def checkout_git_ref(git_ref, verbose=False, repo=None, dry_run=False):
###############################################################################
    """
    Checks out 'branch_ref', and updates submodules
    """
    if dry_run:
        print("Would checkout {}".format(git_ref))
    elif get_current_commit(repo=repo) != get_current_commit(repo=repo, commit=git_ref):
        expect(is_repo_clean(repo=repo), "If we need to change HEAD, then the repo must be clean before running")
        expect(git_ref is not None, "Missing git-ref")

        run_cmd_no_fail("git checkout {}".format(git_ref), from_dir=repo)
        update_submodules(repo=repo)
        git_commit = get_current_commit(repo=repo)
        expect(is_repo_clean(repo=repo), "Something went wrong when checking out git ref '{}'".format(git_ref))

        if verbose:
            print("Switched to '{}' ({})".format(git_ref,git_commit))
            print_last_commit(repo=repo, git_ref=git_ref)

###############################################################################
def get_git_toplevel_dir(repo=None):
###############################################################################
    """
    Get repo toplevel directory. Return None if it could not be found
    """
    stat, output, _ = run_cmd("git rev-parse --show-toplevel", from_dir=repo)
    return output if stat == 0 else None

###############################################################################
def cleanup_repo(orig_branch, orig_commit, repo=None, dry_run=False):
###############################################################################
    """
    Discards all unstaged changes, as well as untracked files
    """
    curr_commit = get_current_commit(repo=repo)

    # Is this a pointless check? Maybe.
    if not dry_run and not is_repo_clean(repo=repo):
        # Discard any modifications to the repo (either tracked or untracked),
        # but keep the ctest-build directory
        run_cmd_no_fail("git clean -df --exclude=ctest-build", from_dir=repo)
        toplevel_dir = get_git_toplevel_dir(repo=repo)
        run_cmd_no_fail("git checkout -- {}".format(toplevel_dir), from_dir=repo)

    checkout_git_ref(orig_branch if orig_branch else orig_commit, repo=repo, dry_run=dry_run)

    # This *can* happen. test_all_scream can merge origin/master into current branch.
    # Checking out orig_branch doesn't do anything if we were on a branch (not detached
    # head mode), since the branch tip moved with the master merge. In that case,
    # what we really need is a hard reset to the original commit.
    # NOTE: if you reset the branch, don't forget to re-update the modules!!
    if curr_commit != orig_commit and not dry_run:
        run_cmd_no_fail("git reset --hard {}".format(orig_commit), from_dir=repo)
        update_submodules(repo=repo)
