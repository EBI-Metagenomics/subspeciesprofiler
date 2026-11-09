# AGENTS.md

<!--
Team rules load from the mi-agentic plugin (Claude) or from ~/.codex/AGENTS.md (Codex).
Do not copy the team rules here. Write only the rules for this repository.
Keep this file short. Delete the sections that you do not need.
-->

## Repository

- Purpose: Subspecies profiling using PopPunk and Syntracker tools
- Main parts: Nextflow pipeline in `workflows/`, tools pipelines in `subworkflows/`, Python package in `bin/`

## Commands

- Run the tests: `nf-test test` and `uv run pytest`
- Run the linters: `nf-core pipelines lint` and `uv run ruff check`

## Rules for this repository

- Never do git commit even in auto mode
- Never remove files from the repo without asking first
