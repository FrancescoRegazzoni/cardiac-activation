# `MATLAB` codes #

## Requirements

- [`MATLAB`](https://www.mathworks.com/products/matlab.html) R2019a or higher

## Usage

Choose the desired model at lines 7-9 of the file [main.m](main.m) and execute it.
This produces a figure with the outcome of the simulation, a print of this figure (named `output_<MODEL_NAME>.pdf`), and a `csv` file named `output_<MODEL_NAME>.csv`, with the following columns:

| column | description                         | measure unit |
|:-------|:------------------------------------|-------------:|
| `t`    | time                                | s            |
| `Ca`   | intracellular calcium concentration | μM           |
| `SL`   | sarcomere length                    | μm           |
| `P`    | permissivity                        | -            |
| `Ta`   | active tension                      | kPa          |