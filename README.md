# About

This repository contains the implementation that performed the experiments for the 
research paper:

* Steven Holtzen, Todd Millstein, and Guy Van den Broeck. ["Generating and
  Sampling Orbits for Lifted Probabilistic
  Inference."](https://arxiv.org/abs/1903.04672) Uncertainty in Artificial
  Intelligence (UAI). 2019.

To cite, please use:

```
@inproceedings{HoltzenUAI19,
  author    = {Holtzen, Steven and Millstein, Todd and Van den Broeck, Guy},
  title     = {Generating and Sampling Orbits for Lifted Probabilistic Inference},
  booktitle = {Proceedings of the 35th Conference on Uncertainty in Artificial Intelligence (UAI)},
  month     = {jul},
  year      = {2019},
}
```

# Disclaimer

This code is provided as-is. It is what was used to run the experiments in the
paper. It is typical research-ware; there are many improvements that could be
made, and it is definitely not ready for prime-time (i.e., production). Think of
it like extra-detailed documentation on how the algorithm described in the paper
works.

# Installation

While it may look like regular Python, this code must actually be run
using the [Sage math library](https://www.sagemath.org/library.html). 
To use this library:

1. Install sage.
2. Install `bliss` by running `sage -i bliss` (if you want to hack the code to
   remove this requirement, it is possible to do so, but it will be much slower)
3. Run `sage setup.py`

# Organization

* `factor.py` contains the most interesting code. It has a simple implementation of 
  a factor graph and the exact and approximate lifted inference algorithms
  from the paper.
* `my_graphs.py` contains the example factor graphs that are generated
* `experiments.py` contains stubs for running the experiments from the paper
* `test.py` contains some standard test cases that illustrate usage.

To execute these files, use `sage` (*not* the regular Python command). I.e.,
run `sage test.py`. 

