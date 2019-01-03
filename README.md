# Linnea

Linnea is an experimental tool for the automatic generation of optimized code for linear algebra problems. It is developed at the [High-Performance and Automatic Computing group](http://hpac.rwth-aachen.de) at RWTH Aachen University.

## Installation

Linnea requires Python 3.6 and can be installed with `pip install git+git://github.com/HPAC/linnea.git`. (Depending on your setup, you may have to use `pip3` instead of `pip`.) To uninstall Linnea, use `pip uninstall linnea`. This also removes the commandline tool.

### Development Installation

If you intend to contribute to Linnea, you can install it from local sources by running `pip install -e .` in your local development directory.

## Overview

Linnea is a prototype of a compiler/program synthesis tool that automates the translation of the mathematical description of a linear algebra problem to an efficient sequence of calls to BLAS and LAPACK kernels. The main idea of Linnea is to construct a search graph that represents a large number of programs, taking into account knowledge about linear algebra, numerical linear algebra and high-performance computing. The algebraic nature of the domain is used to reduce the size of the search graph, without reducing the size of the search space that is explored.

The input to Linnea are linear algebra expressions. As operands, matrices, vectors and scalars are supported. Operands can be annotated with properties, such as 'lower triangular' or 'symmetric'. Supported operations are addition, multiplication, transposition and inversion. At the moment, Linnea generates Julia code (see https://julialang.org), using BLAS and LAPACK wrappers whenever possible.

## Usage

Linnea can be used in two different ways.

### Python Module

At the moment, Linnea is primarily a Python module. An example script for how to use Linnea within Python can found in `examples/run_linnea.py`. The input expressions are represented as Python objects. As an example, consider the description of a lower triangular linear system (omitting imports):

```python
n = 1000

L = Matrix("L", (n, n))
L.set_property(properties.LOWER_TRIANGULAR)
L.set_property(properties.FULL_RANK)
x = Vector("x", (n, 1))
y = Vector("y", (n, 1))

input = Equations(Equal(y, Times(Inverse(L), x)))
```

Further examples of input problems are provided in the `examples/inputX.py` files.

Options can be set with a number of `linnea.config.set_X()` functions.

### Commandline Tool

When installing Linnea via `pip`, the commandline tool `linnea` is installed. As input, it takes a description of the input problem in a simple custom language. With this language, the same lower triangular system is described as:

```
n = 1000

Matrix L(n, n) <LowerTriangular, FullRank>
ColumnVector x(n) <>
ColumnVector y(n) <>

y = inv(L)*x
```

Further examples are provided in `examples/inputX.la`. Notice that the primary purpose of this input format is to make it slightly easier to try out Linnea. There are no plans to establish this as an actual language. New features will probably not be immediately available in this language, and the language may change in the future without being backward compatible.

The list of commandline options is available via `linnea -h`.

### Output

As output, Linnea generates a directory structure that contains code files, as well a file containing a description of the derivation graph, the primary datastructure used by Linnea. Which files are generated can be set as options. Likewise, the location of the output can be specified. By default, it is the current directory.

For the linear system from the previous examples, the following code will be generated:

```julia
using LinearAlgebra.BLAS
using LinearAlgebra

function algorithm0(ml0, ml1)
    # cost 1e+06
    # L: ml0, full, x: ml1, full
    # tmp1 = (L^-1 x)
    trsv!('L', 'N', 'N', ml0, ml1)

    # tmp1: ml1, full
    # y = tmp1
    return (ml1)
end
```

### Options

Linnea offers a number of options which can be set through `linnea.config` in Python or as commandline options for the commandline tool. Alternatively, all options can also be specified in a `linnea_config.json` file (see `examples`) which has to be located in the same directory where Linnea is run, or at the user's `$HOME` folder. Both commandline options and `linnea.config` options override what is specified in `linnea_config.json`. As a fallback, reasonable default values are used.

There are the following options (those are the names used in Python, the commandline options have slightly different names. See `linnea -h`):

* `output_code_path` The output of Linnea will be stored in this directory. The default is the current directory.

* `output_name` Linnea creates a new directory that contains all output files. This is the name of this directory. The default is `tmp`.

* `language` Not available for the commandline tool. For now, the only allowed option is `Julia`.

* `julia_data_type` The data type used in the generated code. Either `Float32` or `Float64`. The default is `Float64`.

* `merging_branches` Whether or not to merge branches in the derivation graph. The default is `true`.

* `dead_ends` Whether or not to eliminate dead ends in the derivation graph early. The default is `true`.

* `algorithms_limit` The upper limit for the number of algorithms that are written to files. The default is `100`.

* `strategy` The strategy used to find algorithms. Either `constructive` or `exhaustive`. The default is `constructive`.

* `generate_graph` Whether or not to generate a `.gv` file of the derivation graph. The default is `false`.

* `graph_style` Style of the derivation graph. Either `full`, `simple`, or `minimal`. The default is `full`. Only applies if `generate_graph` is set to `True`.

* `generate_derivation` Whether or not to generate a description of how the algorithms were derived. The default is `false`.

* `generate_code` Whether or not to generate the actual code of the algorithms. The default is `true`.

* `generate_experiments` Whether or not to generate code that can be used to run the algorithms. The default is `false`.

* `verbosity` Level of verbosity. The default is `1`.

## Publications

A number of publications that discuss different aspects of Linnea can be found [here](http://hpac.rwth-aachen.de/publications/author/Barthels).

## Contributors

* Henrik Barthels
* Marcin Copik
* Diego Fabregat Traver
* Julius Hohnerlein
* Manuel Krebber
