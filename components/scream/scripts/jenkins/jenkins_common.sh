#!/bin/bash -xe

# Check if the env var PR_LABELS is defined, and contains something meaningful
IFS=';' read -r -a labels <<< "$PR_LABELS";

# default values
skip_testing=0
if [ ${#labels[@]} -gt 0 ]; then
  # We do have some labels. Look for some supported ones.
  for label in "${labels[@]}"
  do
    if [ "$label" == "CI: Integrate Without Testing" ]; then
      skip_testing=1
    fi
  done
fi

if [ $skip_testing -eq 0 ]; then
  # User did not request to skip tests. Proceed with testing.
  source ./scream/components/scream/scripts/jenkins/${NODE_NAME}_setup

  git config --global user.email "jenkins@ignore.com"
  git config --global user.name "Jenkins Jenkins"

  SUBMIT="--submit"
  if [ -n "$PULLREQUESTNUM" ]; then
      SUBMIT="" # We don't submit AT runs
  fi

  ./scream/components/scream/scripts/gather-all-data "./scripts/test-all-scream \$compiler -p -i -m \$machine $SUBMIT" -l -m $SCREAM_MACHINE
else
  echo "Tests were skipped, since the Github label 'CI: Integrate Without Testing' was found.\n"
fi
