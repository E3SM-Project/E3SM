(omega-dev-docs)=

# Documentation

The OMEGA documentation is generated using the
[Sphinx](https://www.sphinx-doc.org/en/master/) package and is written in
[MyST](https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html)
format.  We recommend these 
[basic tips](https://myst-parser.readthedocs.io/en/latest/syntax/roles-and-directives.html#roles-directives)
on using MyST in Sphinx.

Another easy way to get started is by taking a look at the existing source
code for the documentation: 

<https://github.com/E3SM-Project/Omega/tree/develop/components/omega/doc>

Each time you add new features to OMEAG, the corresponding documentation must 
be included with the pull request to add the code.  This includes documentation
for both the User's Guide and the Developer's Guide. We will add some examples 
as the documentation gets fleshed out.

Please add anchors for internal linking to the top of every page before
headings within the page whenever you think internal linking might be useful
(now or in the future).  In the User's Guide, please start each anchor with
`omega-`:

```markdown
(omega-quick-start)=

# Quick Start for Users
...
```

This is in anticipation of a future time when the OMEGA documentation might
be combined with other components of E3SM.

In the Developer's Guide, anchors should start with `omega-dev-`:

```markdown
(omega-dev-quick-start)=

# Quick Start for Developers
...
```

Documentation for an OMEGA feature in the User's Guide should contain 
information that is needed for users who set up and run OMEGA, including:

- Flags related to the feature that a user can (or must) set when building
  OMEGA,
- Config options that a user can modify related to the feature

Documentation for the feature in the Developer's Guide may contain:

- Relevant technical details about development specific to that feature
- Details on unit tests that cover that feature
- Details on polaris regression tests that cover that feature
