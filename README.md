# Guido
Guido was developed as an API tool that can search for gRNA targets in any reference genome or DNA sequence. It integrates MMEJ prediction and scoring, off-target search, and allows users to define their own layers of data that can be used in the gRNA evaluation.

## Installation
Todo

## Usage
Todo

## Docs
Todo

## Developer setup
Install [poetry](https://python-poetry.org/docs/#installation):

```bash
$ pip install poetry
```

Create development environment:

```bash
$ cd guido
$ poetry install
```

Activate development environment:

```bash
$ poetry shell
```

Install pre-commit hooks:

```bash
$ pre-commit install
```

Run pre-commit checks (isort, black, blackdoc, flake8, ...) manually:

```bash
$ pre-commit run --all-files
```

Bump version, build and publish to PyPI:

```bash
$ poetry version prerelease
$ poetry build
$ poetry publish
```
