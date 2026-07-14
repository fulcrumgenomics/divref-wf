# Development and Testing

## Primary Development Commands

To check and resolve linting issues in the codebase, run:

```console
uv run --directory divref ruff check --fix
```

To check and resolve formatting issues in the codebase, run:

```console
uv run --directory divref ruff format
```

To check the unit tests in the codebase, run:

```console
uv run --directory divref pytest
```

To check the typing in the codebase, run:

```console
uv run --directory divref mypy
```

To generate a code coverage report after testing locally, run:

```console
uv run --directory divref coverage html
```

To check the `uv` lock file is up to date:

```console
uv lock --directory divref --check
```

To check the `pixi` lock file is up to date:

```console
pixi run check-lock
```

## Shortcut Task Commands

###### For Running Individual Toolkit Checks

```console
uv run --directory divref poe check-format
uv run --directory divref poe check-lint
uv run --directory divref poe check-tests
uv run --directory divref poe check-typing
```

###### For Running All Toolkit Checks

```console
pixi run check-toolkit
```

###### For Running Individual Fixes

```console
uv run --directory divref poe fix-format
uv run --directory divref poe fix-lint
```

###### For Running All Fixes

```console
uv run --directory divref poe fix-all
```

###### For Running All Fixes and Checks

```console
pixi run fix-and-check-all
```
