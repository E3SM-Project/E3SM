from CIME.utils import run_cmd, run_cmd_no_fail, expect, get_timestamp

import getpass, logging, os

# Constants
ESMCI_REMOTE_NAME = "esmci_remote_for_split"
ESMCI_URL = "git@github.com:ESMCI/CIME.git"
SPLIT_TAG_PREFIX = "acme-split-"
MERGE_TAG_PREFIX = "to-acme-"

###############################################################################
def setup():
###############################################################################
    run_cmd_no_fail("git config merge.renameLimit 999999")
    run_cmd_no_fail("git checkout master && git pull", verbose=True)

    remotes = run_cmd_no_fail("git remote")
    if ESMCI_REMOTE_NAME not in remotes:
        run_cmd_no_fail("git remote add {} {}".format(ESMCI_REMOTE_NAME, ESMCI_URL), verbose=True)

    run_cmd_no_fail("git fetch --prune {}".format(ESMCI_REMOTE_NAME), verbose=True)
    run_cmd_no_fail("git fetch --prune {} --tags".format(ESMCI_REMOTE_NAME), verbose=True)

    run_cmd_no_fail("git clean -fd", verbose=True)

###############################################################################
def get_tag(prefix, expected_num=1):
###############################################################################
    tags = run_cmd_no_fail("git tag").split()
    tags = [tag for tag in tags if tag.startswith(prefix)]

    if expected_num == 1:
        return tags[-1]
    else:
        return tags[-expected_num:]

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
    new_tag = "{}{}".format(prefix, get_timestamp(timestamp_format="%Y-%m-%d"))
    expect(old_tag != new_tag, "New tag must have different name than old tag")

    run_cmd_no_fail("git tag {} {}".format(new_tag, commit), verbose=True)
    run_cmd_no_fail("git push {} {}".format(remote, new_tag), verbose=True)

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
def do_subtree_split(new_split_tag, merge_tag):
###############################################################################
    subtree_branch = get_branch_from_tag(new_split_tag)
    run_cmd_no_fail("git subtree split --prefix=cime --onto={} -b {}".\
                        format(merge_tag, subtree_branch), verbose=True)
    return subtree_branch

###############################################################################
def touches_file(start_range, end_range, filepath):
###############################################################################
    return run_cmd_no_fail("git log {}..{} {}".format(start_range, end_range, filepath)) != ""

###############################################################################
def has_non_local_commits(filepath, non_local_path, local_tag):
###############################################################################
    most_recent = get_tag(local_tag)
    return run_cmd_no_fail("git diff MERGE_HEAD:{} {}:{}".format(non_local_path, most_recent, filepath)) != ""

###############################################################################
def reset_file(version, srcpath, dstpath):
###############################################################################
    os.remove(dstpath)
    run_cmd_no_fail("git show {}:{} > {}".format(version, srcpath, dstpath))
    run_cmd_no_fail("git add {}".format(dstpath))

###############################################################################
def get_last_merge(branch_name):
###############################################################################
    return run_cmd_no_fail("git log --first-parent ORIG_HEAD --grep='{}' -1 --oneline".format(branch_name)).split()[0]

###############################################################################
def handle_easy_conflict(filepath, is_merge):
###############################################################################
    non_local_tag = MERGE_TAG_PREFIX if is_merge else SPLIT_TAG_PREFIX
    local_tag     = MERGE_TAG_PREFIX if non_local_tag == SPLIT_TAG_PREFIX else SPLIT_TAG_PREFIX

    local_branch = "branch-for-{}".format(non_local_tag)

    non_local_path = filepath.replace("cime/", "", 1) if is_merge else os.path.join("cime", filepath)

    last_merge = get_last_merge(local_branch)

    if not touches_file(last_merge, "ORIG_HEAD", filepath):
        logging.info("File '{}' appears to have had no recent local mods, setting to non-local".format(filepath))
        reset_file("MERGE_HEAD", non_local_path, filepath)
        return True
    elif not has_non_local_commits(filepath, non_local_path, local_tag):
        logging.info("File '{}' appears to have had no recent non-local mods, setting to local".format(filepath))
        reset_file("ORIG_HEAD", filepath, filepath)
        return True
    else:
        logging.info("File '{}' appears to have real conflicts".format(filepath))
        return False

###############################################################################
def handle_easy_conflicts(is_merge):
###############################################################################
    conflicting_files = run_cmd_no_fail("git diff --name-only --diff-filter=U").splitlines()
    if not conflicting_files:
        expect(False, "Merge appears to have failed for reasons other than merge conflicts")

    rv = []
    for conflicting_file in conflicting_files:
        able_to_handle = handle_easy_conflict(conflicting_file, is_merge)
        if not able_to_handle:
            rv.append(conflicting_file)

    return rv

###############################################################################
def handle_conflicts(is_merge=False, auto_conf=False):
###############################################################################
    logging.info("There are conflicts, analyzing...")
    remaining_conflicts = handle_easy_conflicts(is_merge) if auto_conf else True
    if remaining_conflicts:
        expect(False, "There are merge conflicts. Please fix, commit, and re-run this tool with --resume")
    else:
        logging.info("All conflicts were automatically resovled, continuing")
        run_cmd_no_fail("git commit --no-edit")

###############################################################################
def do_subtree_pull(squash=False, auto_conf=False):
###############################################################################
    stat = run_cmd("git subtree pull {} --prefix=cime {} master".format("--squash" if squash else "", ESMCI_REMOTE_NAME),
                   verbose=True)[0]
    if stat != 0:
        handle_conflicts(is_merge=True, auto_conf=auto_conf)

###############################################################################
def make_pr_branch(branch, branch_head):
###############################################################################
    run_cmd_no_fail("git checkout --no-track -b {} {}".format(branch, branch_head), verbose=True)

    return branch

###############################################################################
def merge_branch(branch, squash=False, auto_conf=False):
###############################################################################
    stat = run_cmd("git merge {} -m 'Merge {branch}' -X rename-threshold=25 {branch}".format("--squash" if squash else "", branch=branch),
                   verbose=True)[0]
    if stat != 0:
        handle_conflicts(auto_conf=auto_conf)

###############################################################################
def delete_tag(tag, remote="origin"):
###############################################################################
    run_cmd_no_fail("git tag -d {}".format(tag), verbose=True)
    run_cmd_no_fail("git push {} :refs/tags/{}".format(remote, tag), verbose=True)

###############################################################################
def e3sm_cime_split(resume, squash=False, auto_conf=False):
###############################################################################
    if not resume:
        setup()

        old_split_tag = get_split_tag()

        try:
            new_split_tag = make_new_split_tag(old_split_tag)

            merge_tag = get_merge_tag()

            pr_branch = do_subtree_split(new_split_tag, merge_tag)

            run_cmd_no_fail("git checkout {}".format(pr_branch), verbose=True)
        except:
            # If unexpected failure happens, delete new split tag
            logging.info("Abandoning split due to unexpected failure")
            delete_tag(new_split_tag)
            raise

        # upstream merge, potential conflicts
        merge_branch("{}/master".format(ESMCI_REMOTE_NAME), squash=squash, auto_conf=auto_conf)

    else:
        old_split_tag, new_split_tag = get_split_tag(expected_num=2)
        logging.info("Resuming split with old tag {} and new tag {}".format(old_split_tag, new_split_tag))
        pr_branch = get_branch_from_tag(new_split_tag)

    run_cmd_no_fail("git push -u {} {}".format(ESMCI_REMOTE_NAME, pr_branch), verbose=True)

###############################################################################
def e3sm_cime_merge(resume, squash=False, auto_conf=False):
###############################################################################
    if not resume:
        setup()

        old_merge_tag = get_merge_tag()

        try:
            new_merge_tag = make_new_merge_tag(old_merge_tag)

            pr_branch = make_pr_branch(get_branch_from_tag(new_merge_tag), "origin/master")
        except:
            logging.info("Abandoning merge due to unexpected failure")
            delete_tag(new_merge_tag, remote=ESMCI_REMOTE_NAME)
            raise

        # potential conflicts
        do_subtree_pull(squash=squash, auto_conf=auto_conf)

    else:
        old_merge_tag, new_merge_tag = get_merge_tag(expected_num=2)
        logging.info("Resuming merge with old tag {} and new tag {}".format(old_merge_tag, new_merge_tag))
        pr_branch = get_branch_from_tag(new_merge_tag)

    run_cmd_no_fail("git push -u origin {}".format(pr_branch), verbose=True)

###############################################################################
def abort_split():
###############################################################################
    new_split_tag = get_split_tag()
    pr_branch = get_branch_from_tag(new_split_tag)
    delete_tag(new_split_tag)
    run_cmd_no_fail("git reset --hard origin/master", verbose=True)
    run_cmd_no_fail("git checkout master", verbose=True)
    run_cmd("git branch -D {}".format(pr_branch), verbose=True)

###############################################################################
def abort_merge():
###############################################################################
    new_merge_tag = get_merge_tag()
    pr_branch = get_branch_from_tag(new_merge_tag)
    delete_tag(new_merge_tag, remote=ESMCI_REMOTE_NAME)
    run_cmd_no_fail("git reset --hard origin/master", verbose=True)
    run_cmd_no_fail("git checkout master", verbose=True)
    run_cmd("git branch -D {}".format(pr_branch), verbose=True)
