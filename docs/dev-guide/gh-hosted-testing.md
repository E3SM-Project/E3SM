# Testing on GitHub-hosted runners

The `e3sm-gh` workflow (`.github/workflows/e3sm-gh-testing.yml`) runs CIME tests with
baseline comparison on GitHub's own runners, which are free for public repositories, rather
than on E3SM's self-hosted `ghci` machines. It builds and runs each test inside an
[E3SM container image](https://github.com/E3SM-Project/containers), with the image's own
CIME machine (`e3sm-container`), the same way the containers repository's `e3sm-case`
workflow does.

Hosted runners start empty every time, so everything a self-hosted machine keeps on disk
between runs has to be stored somewhere else:

| What | Self-hosted `ghci` | GitHub-hosted `e3sm-gh` |
| --- | --- | --- |
| Compilers, MPI, libraries | installed on the machine | the container image |
| Inputdata | machine disk | Actions cache, `e3sm-gh-inputdata-<test>-*` |
| Baselines | machine disk (`BASELINE_ROOT`) | Actions cache, `e3sm-gh-baselines-<image>-<arch>-<test>-*` |

## What runs

To start with, one test:
`SMS_Vmoab_Ln5.ne4pg2_oQU480.WCYCL2010NS`, the coupled case on the MOAB driver, which
puts every active component through the coupler. With no `_P` option, CIME gives this grid
one node, which the image's machine defines as the container's cores: 4 ranks per component
on a hosted runner. It runs in these images:

| Image | Arch | Baselines |
| --- | --- | --- |
| `gnu13-openmpi4` | x86_64, aarch64 | compared |
| `intel2024-openmpi4` | x86_64 | not compared: builds and runs only |

Intel's answers depend on which CPU the hosted runner has: its math libraries pick code paths
by CPU, and GitHub's x86 runners come with different CPU models. The same commit and image
have both matched and differed from intel baselines on different runs, so a comparison
would fail at random. Each job summary shows the runner's CPU model. Building with
`-fimf-arch-consistency=true` (or `-fp-model consistent`) on the `e3sm-container_intel`
machine might make intel reproducible enough to compare, at some cost in speed.

The CUDA images are left out: hosted runners have no GPU, so those tests could only be
built, and there would be nothing to compare.

Tests are listed without machine or compiler. Each job runs `create_test` in the image with
`--machine e3sm-container --compiler <compiler>`, the compiler being set next to each image
in the workflow's matrix. The same test runs on any machine with Docker:

```bash
docker run --rm -v "$PWD:/projects/e3sm/work/E3SM" \
  -v /path/to/inputdata:/projects/e3sm/data/inputdata \
  ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 bash -lc '
    . /projects/e3sm/cime/enable.sh && cd /projects/e3sm/work/E3SM/cime/scripts &&
    ./create_test SMS_Vmoab_Ln5.ne4pg2_oQU480.WCYCL2010NS \
      --machine e3sm-container --compiler gnu --wait'
```

`enable.sh` opts into the image's CIME machine. Until #8926 is merged, the workflow also
downloads the P3 lookup table itself, since nothing else does (#5953).

The workflow runs on PRs (all paths except docs), nightly, and by hand.

## Baselines

### How it works

Baselines are kept in the repository's
[Actions cache](https://docs.github.com/en/actions/reference/workflows-and-actions/dependency-caching),
one entry per image, arch and test. Baselines differ between compilers and between
architectures, so each combination has its own entry.

- **Compare (PRs, nightly).** Each run restores the newest entry for its test and passes
  `-c -b master` to `create_test`. If there is no entry, the test runs without comparing and
  the job shows a warning, so the first run (or one after an entry has expired) does not fail.
- **Bless (by hand, master only).** Go to Actions → `e3sm-gh` → *Run workflow* on `master`
  with *bless* checked. It runs with `-o -g -b master` and saves a new entry, keyed by the
  commit and the run. The next run anywhere restores it. The workflow refuses to bless from
  any branch but master.

This matches how the `ghci` machines are used: baselines change only when someone blesses
on master, such as after merging a PR that is intentionally not bit-for-bit.

### Why the Actions cache

Of the four places considered for baselines, the cache needs the least new machinery:
`actions/cache` with no extra permissions, tokens, or tools. Three of its rules shape the
design:

1. **Entries can't be overwritten.** A bless adds a new entry with a unique key, and every
   run restores the newest one by key prefix.
2. **A PR sees only entries saved on its own ref or on master.** Blessing on master makes
   baselines visible to every PR, and nothing a PR saves can replace them.
3. **Entries unused for 7 days are deleted, and the oldest go first once the repository's
   cache passes its size limit (10 GB by default).** The nightly run restores every
   baseline entry, which keeps them in use. If an entry is deleted anyway, runs stop
   comparing (with a warning) until the next bless.

The other options were:

| Option | Kept for | Drawbacks |
| --- | --- | --- |
| Bless automatically on every merge to master | (cache) | A non-BFB change that merges becomes the baseline without anyone deciding it should |
| Workflow artifacts from the latest bless run | up to 90 days | Needs a lookup for the right run and an `actions: read` token; still expires |
| Container registry (GHCR), pushed with `oras` | until deleted | Needs `packages: write`, a package created and made public by an org admin, and a cleanup policy |

The registry is where baselines should end up (see below). The cache gets the workflow
running first, with the storage confined to two steps.

### Inputdata and baselines are separate entries

Inputdata for the coupled test is about 4 GB and changes only when the test list or E3SM's
input files do. Baselines for these small grids are much smaller and change with each
bless. In one entry, each bless would save the inputdata again and fill the size limit
quickly, so they are kept apart:

- Inputdata has one entry per test, shared by all images and arches. A run restores the
  newest entry and downloads only what it lacks. A run on master that downloaded something
  saves a new entry and, once that entry is listed in the cache, deletes the older ones.
  Runs on any other ref (PRs, or manual runs on a branch) never save inputdata: their
  entries would be visible only to that ref, would still count against the shared limit,
  and would make the cleanup delete master's entry with nothing to replace it.
- Baselines have one entry per image, arch and test, saved only by a bless.

### Known limits

- The images are pulled by their moving tags (for example `gnu13-openmpi4`). A new image
  with a different compiler or library version can change answers. That shows up as a
  baseline failure on every PR until someone blesses again. The job summary records the
  image digest of each run to make that easy to tell apart from a real regression.
- A bless records the E3SM commit and date in `bless_log` next to the baselines (CIME does
  this for E3SM), and a comparison reports it as "Most recent bless" in `TestStatus.log`.
  The workflow also records the image digest in the baseline entry, and warns when a run
  compares against baselines made with a different image.

## Moving baselines to a registry

The plan is to move baselines to a GHCR package (the registry the images already come from)
once the workflow has proven itself. Packages do not expire, can be tagged and rolled back,
and do not count against the cache limit. Only the two baseline steps in the workflow change.
The test step, the inputdata cache, and the triggers stay as they are.

1. **Create the package.** A bless run on master, with `packages: write` added to its
   permissions, pushes the baseline directory with
   [oras](https://oras.land), e.g. as
   `ghcr.io/e3sm-project/e3sm-baselines:<image>-<arch>-<test>-<date>` plus a moving
   `<image>-<arch>-<test>` tag. An org admin then makes the package public and links it to
   this repository, so PRs (including those from forks) can pull it without credentials.
2. **Run both for a while.** Bless runs save to both the cache and the registry. *Restore
   baselines* tries the registry first and falls back to the cache. This checks that the
   registry path works on real PRs without losing the cache.
3. **Switch over.** Remove the cache save and restore for baselines. Add a cleanup policy
   for old dated tags, like the containers repository's `cleanup-ghcr-packages` workflow.
   Optionally pin the image digest next to each baseline tag, so that a new image prompts a
   bless without failing PRs first.
4. **Rolling back** a bad bless becomes moving the `<image>-<arch>-<test>` tag back to an
   earlier dated tag. With the cache, it needs another bless from an earlier commit.

The same approach fits inputdata, which is what limits how many tests the cache can hold
(about 4 GB for the first test, against the 10 GB limit). Inputdata could be published as a
registry image, one layer per set of input files, and mounted read-only into the test
container, so tests would no longer download it or share the cache for it.

## Adding tests

Add the test name to `matrix.test` in the workflow. Inputdata and baseline keys include the
test name, so a new test starts with no entries: its first runs download all its inputdata
and skip comparison until master is blessed. Keep the cache limit in mind: the first test's
inputdata already takes about 4 GB of it.

A test with testmods must name its machine and compiler (for example
`TEST.GRID.COMPSET.e3sm-container_gnu.TESTMOD`), so it goes in the matrix's `include`
list once per compiler rather than in `matrix.test`.
