# `data`

Packaged YAML data read via `importlib.resources` in `read_write.py`. See
render documentation for how to add a new mesh.

## Files

```
data
├── input_files.yaml      # maps each mesh to its input files
└── config_overrides.yaml # automatically applied configuration overrides
```

Both files are validated against `configs/Default.yml` when read.
