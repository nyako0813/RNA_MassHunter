# AGENTS.md

## Workspace

- Work only in `/home/nyako/projects/RNA_MassHunter`.
- Do not access, edit, or copy files from `/mnt/c/Users/nyako/Documents/RNA_MassHunter`.

## Editing constraints

- Do not use `apply_patch` in this repository.
- The sandbox `apply_patch` helper is unavailable for files under `/home/nyako/projects/RNA_MassHunter`.
- Edit files using Python scripts, `cat`/heredoc, `sed`, or other standard Linux file operations instead.
- Do not attempt `apply_patch` first and then fall back.
- Preserve UTF-8 encoding and existing line endings.
- After scripted replacements, verify the replacement count and run `git diff --check`.

## Git safety

- Do not commit unless explicitly instructed.
- Do not commit `config.yaml`, config backups, output, logs, caches, or `.venv`.
- Do not run `git reset --hard`, `git clean`, or destructive restore commands without explicit instruction.
