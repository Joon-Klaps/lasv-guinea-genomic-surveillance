#!/usr/bin/env bash

# Run the `spread` command inside the `spreadgltools-env` conda environment
# without activating the environment in the current shell.
TIME="2024-09-19"

# Run in one line to avoid conda temporary-script parsing issues
conda run -n spreadgltools-env spread --tree "$1" --time "${TIME}" --location "location1,location2"

