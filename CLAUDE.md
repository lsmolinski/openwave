# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OpenWave is an open-source subatomic physics simulator to model the formation of matter and energy from energy wave interactions. The project simulates phenomena from spacetime emergence through particle formation to complex matter behavior.

### Project Goals

To develop OpenWave, an open-source computer simulator with objectives described in the text below, based on the papers attached as 9 files, with special attention to the file `Relationship of the Speed of Light to Aether Density` where there is a Planck mass correction from previous papers (affecting granule mass), using the `06. Constants and Equations - Waves.pdf` as constants reference, built in phases. Simulation physics, constants, and equations will be drawn from research papers located at the `/scientific_source` folder. For performance on the granular physics simulations we'll be using the Taichi Lang python library.

### What is OpenWave?

- Refer to `README.md` for a detailed description and scope of OpenWave.
- Refer to `WELCOME.md` for a quick intro to OpenWave.

### Known Challenges & Limitations

#### Granularity vs. Performance

- Full Planck-scale fidelity may be computationally prohibitive; require user-tunable resolution.
- This project will use a dedicated physics computational backend, independent of 3D modeling software.

## Project Architecture

### Modules Structure and Objects Map

- Refer to `README.md` file for the Modules Structure, Objects Map and System Architecture.

## Installation

- Refer to `README.md` for installation guidance of OpenWave.

## Usage (Work in Progress)

- Refer to `README.md` for usage instructions.

## Scientific Documentation & Requirements

### Scientific Source Material

Each method directory under `openwave/xperiments/` (e.g. `m3_wolff_lafreniere/`, `m5_lagrangian_field/`) contains a `/research` subfolder with on-going research supporting that method.

The `/scientific_source` directory contains foundational papers.

## Physics Context

This project implements Energy Wave Theory concepts:

- Energy Wave as fundamental building blocks
- Wave interactions forming particles and matter
- Simulation from Planck scale to macroscopic phenomena
- Validation against experimental observations

- Also refer to ../CLAUDE.md file to search for any available context to the OpenWave project in a parent directory.

## Code Style & Documentation Standards

- Follow the [Markdown Style Guide](/dev_docs/MARKDOWN_STYLE_GUIDE.md) for all `.md` files
- Adhere to [Coding Standards](/dev_docs/CODING_STANDARDS.md) for Python code
- Apply [Performance Guidelines](/dev_docs/PERFORMANCE_GUIDELINES.md) for optimization
- Use [Loop Optimization](/dev_docs/LOOP_OPTIMIZATION.md) patterns for critical loops

### Important: Markdown Linting Requirements

When editing any `.md` files, ALWAYS ensure compliance with markdown linting rules:

- Add blank lines around headings (before and after)
- Add blank lines around lists (before and after)
- Add blank lines around code blocks (before and after)
- Use consistent ordered list numbering (1, 1, 1 style)
- Check for proper spacing and formatting

Run linting checks after editing to catch any issues before committing.
