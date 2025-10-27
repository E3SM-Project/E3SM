<!--
Please add a description of what is accomplished in the PR here at the top:
-->

<!--
Below are a few things we ask you or your reviewers to kindly check.
***Remove checks that are not relevant by deleting the line(s) below.***
-->
Checklist
* [ ] Documentation:
  * [ ] Design document has been generated and [added to the docs](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/design)
  * [ ] [User's Guide](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/userGuide) has been updated
  * [ ] [Developer's Guide](https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc/devGuide) has been updated
  * [ ] Documentation has been [built locally](https://docs.e3sm.org/Omega/omega/develop/devGuide/BuildDocs.html) and changes look as expected
* [ ] Linting
  * [ ] [Pre-commit](https://docs.e3sm.org/Omega/omega/develop/devGuide/Linting.html) run locally
* [ ] Building
  * [ ] CMake build does not produce any new warnings from changes in this PR
* [ ] Testing
  * [ ] Add a comment to the PR titled `Testing` with the following:
    * [ ] Which machines [CTest unit tests](https://docs.e3sm.org/Omega/omega/develop/devGuide/QuickStart.html#running-ctests)
          have been run on and indicate that are all passing.
    * [ ] The [Polaris omega_pr test suite](https://docs.e3sm.org/Omega/omega/develop/devGuide/Testing.html)
          has passed, using `develop` as a baseline
    * [ ] Document machine(s), compiler(s) and test path(s) (`-p` in `polaris suite`) used
    * [ ] Indicate "All tests passed" or document failing tests
    * [ ] Document testing used to verify the changes including any tests that are added/modified/impacted.
    * [ ] Performance related PRs: Please include a relevant PACE experiment link documenting performance before and after.
  * [ ] New tests:
    * [ ] CTest unit tests for new features have been added per the approved design.
    * [ ] Polaris tests for new features have been added per the approved design (and included in a test suite)
* [ ] Stealth Features
  * [ ] If any stealth features are included in the PR, please confirm that they have been documented.

<!--
Please note any issues this fixes using closing keywords: https://help.github.com/articles/closing-issues-using-keywords
-->


