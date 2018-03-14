# Contributing

## Introduction
First off, thank you for considering contributing to CIME. CIME is a community-driven
project, so it's people like you that make CIME useful and successful.

Following these guidelines helps to communicate that you respect the time of the
developers managing and developing this open source project. In return, they
should reciprocate that respect in addressing your issue, assessing changes, and
helping you finalize your pull requests.

We love contributions from community members, just like you! There are many ways
to contribute, from writing tutorials or examples, improvements
to the documentation, submitting bug report and feature requests, or even writing
code which can be incorporated into CIME for everyone to use. If you get stuck at
any point you can create an [issue on GitHub](https://github.com/ESMCI/CIME/issues).

For more information on contributing to open source projects,
[GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point. Also, checkout the [Zen of Scientific Software Maintenance](https://jrleeman.github.io/ScientificSoftwareMaintenance/)
for some guiding principles on how to create high quality scientific software contributions.

## Getting Started

Interested in helping extend CIME? Have code from your research that you believe others will
find useful? Have a few minutes to tackle an issue? In this guide we will get you setup and
integrated into contributing to CIME!

## What Can I Do?
* Tackle any unassigned [issues](https://github.com/ESMCI/CIME/issues) you wish!  

* Contribute code you already have. It doesn’t need to be perfect! We will help you clean
  things up, test it, etc.

* Make a tutorial or example of how to do something.

## How Can I Talk to You?
Discussion of CIME development often happens in the issue tracker and in pull requests.

## Ground Rules
The goal is to maintain a diverse community that's pleasant for everyone.

* Each pull request should consist of a logical collection of changes. You can
  include multiple bug fixes in a single pull request, but they should be related.
  For unrelated changes, please submit multiple pull requests.
* Do not commit changes to files that are irrelevant to your feature or bugfix
  (eg: .gitignore).
* Be willing to accept constructive criticism as part of issuing a pull request,
  since the CIME developers are dedicated to ensuring that new features extend the
  system robustly and do not introduce new bugs.
* Be aware that the pull request review process is not immediate, and is
  generally proportional to the size of the pull request.

## Reporting a bug
When creating a new issue, please be as specific as possible. Include the version
of the code you were using, as well as what operating system you are running.
If possible, include complete, minimal example code that reproduces the problem.

## Pull Requests
**Working on your first Pull Request?** You can learn how from this *free* video series [How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github) or the guide [“How to Contribute to Open Source"](https://opensource.guide/how-to-contribute/).
We love pull requests from everyone. Fork, then clone the repo:

    git clone git@github.com:your-username/CIME.git

You will need to install CIME dependencies and edit config files
to tell CIME about your development machine. See the [CIME users guide](https://esmci.github.io/cime/users_guide/porting-cime.html)

Run the scripts_regression_test:

    cd scripts/tests
    scripts_regression_tests.py

Make your change. Add tests for your change. Make the tests pass to the same level as before your changes.

    scripts_regression_tests.py

Commit the changes you made. Chris Beams has written a [guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages.

Push to your fork and [submit a pull request][pr].

[pr]: https://github.com/ESMCI/CIME/compare

## Code Review
Once you have submitted a pull request, expect to hear at least a comment within a couple of days.
We may suggest some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted:

* Write tests.
* Follow [PEP8][pep8] for style. (The `flake8` utility can help with this.)
* Write a [good commit message][commit].

Pull requests will automatically have tests run by Travis. This includes
running both the unit tests as well as the `flake8` code linter.

[pep8]: http://pep8.org
[commit]: https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
