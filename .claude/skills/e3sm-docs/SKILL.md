---
name: e3sm-docs
description: Use when creating, editing, reviewing, or planning E3SM Markdown/MkDocs documentation under docs/, components/*/docs/, tools/*/docs/, or related mkdocs.yml navigation files.
---

# E3SM documentation authoring skill (Claude-compatible copy)

## Purpose

When invoked, provide concise, repo-aware guidance for creating, editing, and reviewing
Markdown documentation for E3SM. Use this skill for docs under `docs/`,
`components/*/docs/`, and `tools/*/docs/`, and for changes to `mkdocs.yaml` or
component `mkdocs.yml` files.

## Maintenance requirement

This skill is intentionally duplicated in two repository locations:

- `.agents/skills/e3sm-docs/SKILL.md`
- `.claude/skills/e3sm-docs/SKILL.md`

When editing either copy, make the equivalent change to the other copy in the SAME commit or PR.
Before finalizing the change, compare both files and confirm they are identical, unless there is an
explicitly documented, tool-specific reason for divergence.

If you are an automated tool or bot editing these files, ensure your change updates both paths and
include a short note in the commit message: "sync: update e3sm-docs skill (agents + claude)".

## Repo documentation layout

- Top-level site: `docs/`
- Component docs: `components/<component>/docs/`
- Tool docs: `tools/<tool>/docs/`
- Site config: `mkdocs.yaml` (root) and `components/*/mkdocs.yml` (component-specific nav)

## Authoring rules (concise)

- Use relative links for repo-local references.
- Use fenced code blocks with language tags for command and code examples.
- Use MkDocs Material admonitions (!!! note/warning) where helpful.
- Keep changes factual and cite source files (configs, scripts, XML) when asserting behavior.
- Do not invent scientific claims, defaults, or validation results. Ask for domain expert review when unsure.
- Minimize broad markdownlint disables; prefer inline, single-rule disables with justification.

## Validation checklist

1. Run markdownlint if available (local or CI) using `docs/.markdownlint.json` for rules.
2. For navigation changes, update the appropriate `mkdocs.yml` and ensure links are correct.
3. CI runs `mkdocs build --strict --verbose`; note that some dev machines may not have mkdocs installed.
4. Add or update the PR description to mention docs changes and any expected CI failures.

## Quick example: adding a new docs page

1. Create `components/<component>/docs/<new>.md` with required frontmatter if needed.
2. Update the component `mkdocs.yml` to include the new page in navigation.
3. Run markdownlint (locally if possible) and note any rules you intentionally relax.
4. Open a PR that updates both skill copies if you modified `SKILL.md`.

## Contacts

When in doubt about scientific intent or intended model behavior, ask the component maintainers or
file an issue and @-mention the relevant team.
