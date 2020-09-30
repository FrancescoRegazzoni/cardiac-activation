# `Python` codes #

## Requirements

- [`Python`](https://www.python.org/) 3.7 or higher
- [`NumPy`](https://numpy.org/) 1.18 or higher
- [`SciPy`](https://www.scipy.org/) 1.5 or higher
- [`Matplotlib`](https://matplotlib.org/) 3.1.3 or higher (only used to plot the results)
- [`pandas`](https://pandas.pydata.org/) 1.0.0 or higher (only used to save the results into a `csv` file)

## Usage

Choose the desired model at lines 15-17 of the file [main.py](main.py) and execute it.
This produces a figure with the outcome of the simulation, a print of this figure (named `output_<MODEL_NAME>.pdf`), and a `csv` file named `output_<MODEL_NAME>.csv`, with the following columns:

| column | description                         | measure unit |
|:-------|:------------------------------------|-------------:|
| `t`    | time                                | s            |
| `Ca`   | intracellular calcium concentration | μM           |
| `SL`   | sarcomere length                    | μm           |
| `P`    | permissivity                        | -            |
| `Ta`   | active tension                      | kPa          |