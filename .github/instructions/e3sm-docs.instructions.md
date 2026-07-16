---
applyTo: "docs/**/*.md,components/*/docs/**/*.md,tools/*/docs/**/*.md,**/mkdocs.y*ml,.agents/skills/**/SKILL.md,.claude/skills/**/SKILL.md"
---

When editing E3SM documentation or the repository's documentation skills, keep the following
rules in mind:

- The documentation skill is duplicated at `.agents/skills/e3sm-docs/SKILL.md` and
  `.claude/skills/e3sm-docs/SKILL.md`. Update both files in the same commit or PR and keep them
  identical unless a documented, tool-specific reason exists.
- Use relative links for repo-local references and update the appropriate `mkdocs.yml` when
  adding public pages.
- Run markdownlint (see `docs/.markdownlint.json`) and prefer narrow, justified disables.
- Document validation steps in the PR description. CI runs `mkdocs build --strict --verbose` and
  will report site errors.
