# `C++` codes #

## Requirements

- a `C++` compiler (e.g. [`GCC`](https://gcc.gnu.org/) 7.5 or higher)
- [`CMake`](https://cmake.org/) 3.10 or higher
- [`Boost`](https://www.boost.org/) 1.41 or higher
- [`Eigen`](http://eigen.tuxfamily.org/) 3.3 or higher

## Installation

1. Run the configuration:
```bash
mkdir build
cd build
cmake ..
```

2. Run the compilation:
```bash
make
```

## Usage

To run the code, execute `run_model` passing as unique argument the model name (either `RDQ18`, `RDQ20-MF` or `RDQ20-SE`).
For example:
```bash
./run_model RDQ20-MF
```
This produces a `csv` file named `output_<MODEL_NAME>.csv`, with the following columns:

| column   | description                             | measure unit |
|:---------|:----------------------------------------|-------------:|
| `t`      | time                                    | s            |
| `Ca`     | intracellular calcium concentration     | μM           |
| `SL`     | sarcomere length                        | μm           |
| `dSL_dt` | time derivative of the sarcomere length | μm/s         |
| `Ta`     | active tension                          | kPa          |