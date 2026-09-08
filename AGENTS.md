# nf-cmgg/preprocessing

This repository uses the nf-core pipeline template (`is_nfcore: false` in `.nf-core.yml`) in the nf-cmgg organisation.

## nf-core instructions

Follow [nf-core pipeline AGENTS.md](https://github.com/nf-core/agents/blob/main/resources/pipeline/AGENTS.md). Those rules apply here. Do not restate them in this file.

## This repository

- Use `pixi install` then `pixi shell` for development. `pixi.toml` provides nextflow, nf-core, nf-test, and prek.
- Institutional configs come from [nf-cmgg/configs](https://github.com/nf-cmgg/configs) via `params.custom_config_base` in `nextflow.config`, not from nf-core/configs.
- Main analysis lives in `workflows/preprocessing.nf`. Pipeline-only code belongs in `modules/local/` and `subworkflows/local/`.
- Use `nf-core modules install` to add new modules. Remember to also amend citations across the repository.
- Use `nf-metro` to render custom workflow diagram based on the mermaid source in `docs/images/metro_map_light.md` and `docs/images/metro_map_dark.md`.
- Do not invent files that nf-core tools should generate.

## After writing

After any code or docs edit, run the deslop skill on the new diff. Strip extra comments, defensive noise, and style that does not match the surrounding files. Keep behaviour unchanged.

After Nextflow code edits, run on the files you changed:

```
nextflow lint -format -sort-declarations -exclude ".nf-test" -exclude ".pixi" -exclude "modules/nf-core" -exclude "subworkflows/nf-core" -harshil-alignment .
```

Resolve reported issues. Do not format `modules/nf-core` or `subworkflows/nf-core`.
