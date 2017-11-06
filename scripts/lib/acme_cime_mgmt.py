from CIME.utils import run_cmd, run_cmd_no_fail, expect, get_timestamp

import sys, getpass, logging

# Constants
ESMCI_REMOTE_NAME = "esmci_remote_for_split"
ESMCI_URL = "git@github.com:ESMCI/CIME.git"
SPLIT_TAG_PREFIX = "acme-split-"
MERGE_TAG_PREFIX = "to-acme-"

###############################################################################
def setup():
###############################################################################
    run_cmd_no_fail("git config merge.renameLimit 999999")
    run_cmd_no_fail("git checkout master && git pull && git submodule update --init")

    remotes = run_cmd_no_fail("git remote")
    if ESMCI_REMOTE_NAME not in remotes:
        run_cmd_no_fail("git remote add {} {}".format(ESMCI_REMOTE_NAME, ESMCI_URL))

    run_cmd_no_fail("git fetch {}".format(ESMCI_REMOTE_NAME))
    run_cmd_no_fail("git fetch {} --tags".format(ESMCI_REMOTE_NAME))

###############################################################################
def get_tag(prefix, expected_num=1):
###############################################################################
    tags = run_cmd_no_fail("git tag").split()
    tags = [tag for tag in tags if tag.startswith(prefix)]

    expect(len(tags) == expected_num, "Expected exactly {} {} tag, found {}".format(expected_num, prefix, ", ".join(tags)))

    if expected_num == 1:
        return tags[0]
    else:
        return tags

###############################################################################
def get_split_tag(expected_num=1):
###############################################################################
    return get_tag(SPLIT_TAG_PREFIX, expected_num=expected_num)

###############################################################################
def get_merge_tag(expected_num=1):
###############################################################################
    return get_tag(MERGE_TAG_PREFIX, expected_num=expected_num)

###############################################################################
def make_new_tag(prefix, old_tag, remote="origin", commit="HEAD"):
###############################################################################
    new_tag = "{}{}".format(prefix, get_timestamp(timestamp_format="%m-%d-%Y"))
    expect(old_tag != new_tag, "New tag must have different name than old tag")

    run_cmd_no_fail("git tag {} {}".format(new_tag, commit))
    run_cmd_no_fail("git push {} {}".format(remote, new_tag))

    return new_tag

###############################################################################
def make_new_split_tag(old_split_tag):
###############################################################################
    return make_new_tag(SPLIT_TAG_PREFIX, old_split_tag)

###############################################################################
def make_new_merge_tag(old_merge_tag):
###############################################################################
    return make_new_tag(MERGE_TAG_PREFIX, old_merge_tag,
                        remote=ESMCI_REMOTE_NAME, commit="{}/master".format(ESMCI_REMOTE_NAME))

###############################################################################
def get_branch_from_tag(tag):
###############################################################################
    branch = "{}/branch-for-{}".format(getpass.getuser(), tag)
    return branch

###############################################################################
def do_subtree_split(old_split_tag, new_split_tag, merge_tag):
###############################################################################
    subtree_branch = get_branch_from_tag(new_split_tag)
    run_cmd_no_fail("git subtree split {}.. --prefix=cime --onto={} --ignore-joins -b {}".\
                        format(old_split_tag, merge_tag, subtree_branch))
    return subtree_branch

###############################################################################
def do_subtree_pull():
###############################################################################
    stat = run_cmd("git subtree pull --prefix=cime {} master".format(ESMCI_REMOTE_NAME))[0]
    if stat != 0:
        logging.info("There are merge conflicts. Please fix, commit, and re-run this tool with --resume")
        sys.exit(1)

###############################################################################
def make_pr_branch(branch, branch_head):
###############################################################################
    pr_branch = "{}-pr".format(branch)
    run_cmd_no_fail("git checkout -b {} {}".format(pr_branch, branch_head))

    return pr_branch

###############################################################################
def merge_branch(branch, resume_count):
###############################################################################
    stat = run_cmd("git merge -m 'Merge {}' -X rename-threshold=25 {}".format(branch, branch))[0]
    if stat != 0:
        logging.info("There are merge conflicts. Please fix, commit, and re-run this tool with --resume-{}".format(resume_count))
        sys.exit(1)

###############################################################################
def merge_pr_branch_1(subtree_branch):
###############################################################################
    merge_branch(subtree_branch, "one")

###############################################################################
def merge_pr_branch_2():
###############################################################################
    merge_branch("{}/master".format(ESMCI_REMOTE_NAME), "two")

###############################################################################
def delete_tag(tag, remote="origin"):
###############################################################################
    run_cmd_no_fail("git tag -d {}".format(tag))
    run_cmd_no_fail("git push {} :refs/tags/{}".format(remote, tag))

###############################################################################
def acme_cime_split(resume_one, resume_two):
###############################################################################
    if not resume_one and not resume_two:
        setup()

        old_split_tag = get_split_tag()

        try:
            new_split_tag = make_new_split_tag(old_split_tag)

            merge_tag = get_merge_tag()

            subtree_branch = do_subtree_split(old_split_tag, new_split_tag, merge_tag)

            pr_branch = make_pr_branch(subtree_branch, merge_tag)
        except:
            # If unexpected failure happens, delete new split tag
            delete_tag(new_split_tag)
            raise

        merge_pr_branch_1(subtree_branch)
    else:
        old_split_tag, new_split_tag = get_split_tag(expected_num=2)
        pr_branch = "{}-pr".format(get_branch_from_tag(new_split_tag))

    if not resume_two:
        merge_pr_branch_2()

    try:
        run_cmd_no_fail("git push {} {}".format(ESMCI_REMOTE_NAME, pr_branch))
    except:
        delete_tag(old_split_tag)
        raise

    delete_tag(old_split_tag)

###############################################################################
def acme_cime_merge(resume):
###############################################################################
    if not resume:
        setup()

        old_merge_tag = get_merge_tag()

        try:
            new_merge_tag = make_new_merge_tag(old_merge_tag)

            pr_branch = make_pr_branch(get_branch_from_tag(new_merge_tag), "origin/master")
        except:
            delete_tag(new_merge_tag, remote=ESMCI_REMOTE_NAME)
            raise

        do_subtree_pull()

    else:
        old_merge_tag, new_merge_tag = get_merge_tag(expected_num=2)
        pr_branch = "{}-pr".format(get_branch_from_tag(new_merge_tag))

    try:
        run_cmd_no_fail("git push origin {}".format(pr_branch))
    except:
        delete_tag(old_merge_tag, remote=ESMCI_REMOTE_NAME)
        raise

    delete_tag(old_merge_tag, remote=ESMCI_REMOTE_NAME)
