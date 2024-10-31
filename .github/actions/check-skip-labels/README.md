# Composite action to check if we can skip a job for a PR

This action is meant to be used inside a PR testing workflow, as

```yaml
jobs:
  my-testing:
    steps:
      ...
      - name: check-pr-labels
        if: github.event_name == "pull_request" || github.event_name == "pull_request_review"
        uses: ./.github/actions/check-skip-labels
        with:
          skip_labels: label1,label2,label3
```

The input skip_label is a comma-separated list of labels that, if found
on the PR, will cause this job to terminate immediately with a PASS state.

Ideally, we would like to run this check at the job level, so that we can
skip the job altogether (without using runner time). But while for the
pull_request event we DO have access to the labels from the gh context
(and therefore can check), for pull_request_review we don't, so we need
to ping github for some information
