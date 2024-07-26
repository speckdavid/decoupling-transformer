# Decoupled Search for the Masses: A Novel Task Transformation for Classical Planning

This repository contains the code for the ICAPS 2024 paper ([[pdf]](https://mrlab.ai/papers/speck-gnad-icaps2024.pdf) [[bib]](https://mrlab.ai/papers/speck-gnad-icaps2024.html)). 
It extends [Fast Downward 23.06](https://github.com/aibasel/downward/) by the possibility of running the search on a decoupled representation of the planning problem.

## Build instructions

See [BUILD.md](BUILD.md). For best performance of decoupled search, it is highly recommended to install CPLEX as described in the file. Although a preliminary version without CPLEX is available, CPLEX significantly improves results, so its installation is strongly advised.

## Recommended Search Configurations

We recommend the following search configuration which we defined as aliases.

```
./fast-downward.py --alias "decoupled-lama-first" domain.pddl problem.pddl
```

See [driver/aliases.py](driver/aliases.py) for more aliases.

## Other Search Configuration

You can also create your own search configuration on the decoupled task as follows.

```
./fast-downward.py --root-task-transform "decoupled(factoring=lp())" --search "XXX"
```

Here, the search "XXX" can be selected as in normal Fast Downward, e.g., `"lazy_greedy([ff()])"`. For more information, see the [Fast Downward website](https://www.fast-downward.org).

### Parameters

```
--root-task-transform decoupled(factoring=lp(), same_leaf_preconditons_single_variable=true, conclusive_leaf_encoding=multivalued, skip_unnecessary_leaf_effects=true, conclusive_operators=true, normalize_task=false, normalize_variable_names=false, dump_task=false, write_sas=false, write_pddl=false, write_factoring=false)
```

- `factoring` (Factoring): method that computes the factoring. (see below)
- `same_leaf_preconditons_single_variable` (bool): The same preconditions of leaves have a single secondary variable.
- `conclusive_leaf_encoding` ({basic, binary, multivalued}): Conclusive leaf encoding.
  - `basic`: no special treatment for conclusive leaves. Operators have conditional effects regarding conclusive leaves.
  - `binary`: primary conclusive leaf variables are represented by binary variables. Operators do not have conditional effects regarding a conclusive leaf; instead, they set the primary variable corresponding to the unique reached leaf state to true and all others to false.
  - `multivalued`: primary conclusive leaf variables are represented using the original variables in a factored manner. Operators do not have conditional effects regarding a conclusive leaf; they simply set the primary leaf variables to the corresponding values of the reached leaf state.
- `skip_unnecessary_leaf_effects` (bool): Skip unnecessary leaf effects for operators that have no influence on the leaf.
- `conclusive_operators` (bool): Avoid conditional effects for the effects of conclusive operators on a non-conclusive leaf.
- `normalize_task` (bool): Sort conditions and effects according to variable ids.
- `normalize_variable_names` (bool): Normalizes the variable names by numbering in the format var[x]
- `dump_task` (bool): Dumps the task to the console.
- `write_sas` (bool): Writes the decoupled task to dec_output.sas.
- `write_pddl` (bool): Writes the decoupled task to dec_domain.pddl and dec_problem.pddl.
- `write_factoring` (bool): Writes the factoring of the decoupled task to factoring.txt.

#### LP Factoring

```
lp(verbosity=normal, min_number_leaves=2, max_leaf_size=1000000, factoring_time_limit=30, prune_fork_leaf_state_spaces=false, strategy=MFA, min_mobility=1, min_flexibility=0.2, min_fact_flexibility=0, add_cg_sccs=true)
```

- `verbosity` ({silent, normal, verbose, debug}): Option to specify the verbosity level.
  - `silent`: only the most basic output
  - `normal`: relevant information to monitor progress
  - `verbose`: full output
  - `debug`: like verbose with additional debug output
- `min_number_leaves` (int): The minimum number of leaf factors.
- `max_leaf_size` (int): Maximum domain-size product of variables in a leaf.
- `factoring_time_limit` (int): Time limit for computing the factoring.
- `prune_fork_leaf_state_spaces` (bool): Run simulation-based pruning in fork leaves to reduce their state space, not supported, yet.
- `strategy` ({mml, mmas, mm_opt, mm_approx, mfa, mm, mcl, mcm}): This option determines the property of the factoring that is being optimized by the LP, e.g. the number of mobile leaves, or the sumof leaf mobility.
  - `MML`: maximize mobile leaves
  - `MMAS`: maximize mobile action schemas
  - `MM_OPT`: maximize mobility
  - `MM_APPROX`: maximize mobility (approximation)
  - `MFA`: maximize mobile facts
  - `MM`: maximize mobility (sum)
  - `MCL`: maximize number of mobile conclusive leaves
  - `MCM`: maximize conclusive mobility, i.e. number of conclusive actions (sum)
- `min_mobility` (int): Minimum number of leaf-only actions per leaf factor.
- `min_flexibility` (double): Minimum flexibility (ratio between the number of leaf-only vs. all actions affecting a leaf.
- `min_fact_flexibility` (double): Fact flexibility is measured as the mean ratio across all facts in a leaf of the number ofleaf-only vs. all actions affecting that leaf. This option imposes a minimum on that metric.
- `add_cg_sccs` (bool): If true, every SCC of the causal graph is considered a leaf candidate.

#### Maximum-weight independent set factoring

MWIS factoring is an alternative with a preliminary implementation that does not rely on CPLEX, but it often yields inferior performance. Thus, we recommend to use the LP factoring method. 

```
wmis(verbosity=normal, min_number_leaves=2, max_leaf_size=1000000, factoring_time_limit=30, prune_fork_leaf_state_spaces=false, strategy=MML, min_mobility=1, min_flexibility=0.2, min_fact_flexibility=0, add_cg_sccs=true)
```

- `verbosity` ({silent, normal, verbose, debug}): Option to specify the verbosity level.
  - `silent`: only the most basic output
  - `normal`: relevant information to monitor progress
  - `verbose`: full output
  - `debug`: like verbose with additional debug output
- `min_number_leaves` (int): The minimum number of leaf factors.
- `max_leaf_size` (int): Maximum domain-size product of variables in a leaf.
- `factoring_time_limit` (int): Time limit for computing the factoring.
- `prune_fork_leaf_state_spaces` (bool): Run simulation-based pruning in fork leaves to reduce their state space, not supported, yet.
- `strategy` ({mml, mmas, mm_opt, mfa, mm, mcl, mcm}): This option determines the property of the factoring that is being optimized by the LP, e.g. the number of mobile leaves, or the sumof leaf mobility.
  - `MML`: maximize mobile leaves
  - `MMAS`: maximize mobile action schemas
  - `MM_OPT`: maximize mobility
  - `MFA`: maximize mobile facts
  - `MM`: maximize mobility (sum)
  - `MCL`: maximize number of mobile conclusive leaves
  - `MCM`: maximize conclusive mobility, i.e. number of conclusive actions (sum)
- `min_mobility` (int): Minimum number of leaf-only actions per leaf factor.
- `min_flexibility` (double): Minimum flexibility (ratio between the number of leaf-only vs. all actions affecting a leaf.
- `min_fact_flexibility` (double): Fact flexibility is measured as the mean ratio across all facts in a leaf of the number ofleaf-only vs. all actions affecting that leaf. This option imposes a minimum on that metric.
- `add_cg_sccs` (bool): If true, every SCC of the causal graph is considered a leaf candidate.

#### Miura & Fukunaga factoring

Implementation of Miura and Fukunaga's factoring. It finds a single leaf and is mainly used for comparison. Generally, the LP factoring method supersedes this variant.

```
mf(verbosity=normal, min_number_leaves=2, max_leaf_size=1000000, factoring_time_limit=30, prune_fork_leaf_state_spaces=false)
```

- `verbosity` ({silent, normal, verbose, debug}): Option to specify the verbosity level.
  - `silent`: only the most basic output
  - `normal`: relevant information to monitor progress
  - `verbose`: full output
  - `debug`: like verbose with additional debug output
- `min_number_leaves` (int): The minimum number of leaf factors.
- `max_leaf_size` (int): Maximum domain-size product of variables in a leaf.
- `factoring_time_limit` (int): Time limit for computing the factoring.
- `prune_fork_leaf_state_spaces` (bool): Run simulation-based pruning in fork leaves to reduce their state space, not supported, yet

## Decoupled task to SAS or PDDL

It is possible to use our tool as both a pre-processor and a post-processor to decouple a planning problem, write it to PDDL files or the Fast Downward-specific SAS file, and then transform a decoupled plan into an original one. The following three steps are required.

### Write the decoupled task to a file

The following will write three files: `original.sas`, which represents the original problem; `factoring.txt`, which specifies the factoring; and `dec_output.sas`, which represents the decoupled problem.

```
./fast-downward.py --sas-file original.sas --keep-sas-file domain.pddl problem.pddl --root-task-transform "decoupled(factoring=lp(),write_factoring=true,write_sas=true)" --search "astar(blind(),bound=0)"
```

Or we can write the problem in `dec_domain.pddl` and `dec_problem.pddl`.

```
./fast-downward.py --sas-file original.sas --keep-sas-file domain.pddl problem.pddl --root-task-transform "decoupled(factoring=lp(),write_factoring=true,write_pddl=true)" --search "astar(blind(),bound=0)"
```

### Solve the decoupled task with your own solver

You can now solve the problem using your own solver. Below is an example using Fast Downward on the written decoupled task in `dec_output.sas`.

```
./fast-downward.py --keep-sas-file --plan-file decoupled_plan --alias "lama-first" dec_output.sas
```

### Post-process your plan

Finally, we have a decoupled plan written to `decoupled_plan` which we need to post-process by adding leave actions. This can be done as follows. Note that it is important to input the original `original.sas`.

```
./fast-downward.py original.sas --root-task-transform "decoupled_plan_reconstruction()" --search "astar(blind(),bound=0)"
```

## Citation
David Speck and Daniel Gnad:
Decoupled Search for the Masses: A Novel Task Transformation for Classical Planning.
ICAPS 2024.
[[pdf]](https://mrlab.ai/papers/speck-gnad-icaps2024.pdf) [[bib]](https://mrlab.ai/papers/speck-gnad-icaps2024.html)

---


<img src="misc/images/fast-downward.svg" width="800" alt="Fast Downward">

Fast Downward is a domain-independent classical planning system.

Copyright 2003-2023 Fast Downward contributors (see below).

For further information:
- Fast Downward website: <https://www.fast-downward.org>
- Report a bug or file an issue: <https://issues.fast-downward.org>
- Fast Downward mailing list: <https://groups.google.com/forum/#!forum/fast-downward>
- Fast Downward main repository: <https://github.com/aibasel/downward>

## Scientific experiments

We recommend to use the [latest release](https://github.com/aibasel/downward/releases/latest) instead of the tip of the main branch.
The [Downward Lab](https://lab.readthedocs.io/en/stable/) Python package helps running Fast Downward experiments.
Our separate [benchmark repository](https://github.com/aibasel/downward-benchmarks) contains a collection of planning tasks.

## Supported software versions

The planner is mainly developed under Linux; and all of its features should work with no restrictions under this platform.
The planner should compile and run correctly on macOS, but we cannot guarantee that it works as well as under Linux.
The same comment applies for Windows, where additionally some diagnostic features (e.g., reporting peak memory usage when the planner is terminated by a signal) are not supported.
Setting time and memory limits and running portfolios is not supported under Windows either.

This version of Fast Downward has been tested with the following software versions:

| OS           | Python | C++ compiler                                                     | CMake |
| ------------ | ------ | ---------------------------------------------------------------- | ----- |
| Ubuntu 22.04 | 3.10   | GCC 11, GCC 12, Clang 14                                         | 3.22  |
| Ubuntu 20.04 | 3.8    | GCC 10, Clang 12                                                 | 3.16  |
| macOS 12     | 3.10   | AppleClang 14                                                    | 3.24  |
| macOS 11     | 3.8    | AppleClang 13                                                    | 3.24  |
| Windows 10   | 3.8    | Visual Studio Enterprise 2019 (MSVC 19.29) and 2022 (MSVC 19.31) | 3.22  |

We test LP support with CPLEX 22.1.1 and SoPlex 6.0.3+. On Ubuntu we
test both CPLEX and SoPlex. On Windows we currently only test CPLEX,
and on macOS we do not test LP solvers (yet).

## Build instructions

See [BUILD.md](BUILD.md).


## Contributors

The following list includes all people that actively contributed to
Fast Downward, i.e., all people that appear in some commits in Fast
Downward's history (see below for a history on how Fast Downward
emerged) or people that influenced the development of such commits.
Currently, this list is sorted by the last year the person has been
active, and in case of ties, by the earliest year the person started
contributing, and finally by last name.

- 2003-2023 Malte Helmert
- 2008-2016, 2018-2023 Gabriele Roeger
- 2010-2023 Jendrik Seipp
- 2010-2011, 2013-2023 Silvan Sievers
- 2012-2023 Florian Pommerening
- 2013, 2015-2023 Salomé Eriksson
- 2015, 2021-2023 Thomas Keller
- 2018-2023 Patrick Ferber
- 2018-2020, 2023 Augusto B. Corrêa
- 2021-2023 Clemens Büchner
- 2022-2023 Remo Christen
- 2023 Simon Dold
- 2023 Claudia S. Grundke
- 2023 Emanuele Tirendi
- 2021-2022 Dominik Drexler
- 2016-2020 Cedric Geissmann
- 2017-2020 Guillem Francès
- 2020 Rik de Graaff
- 2015-2019 Manuel Heusner
- 2017 Daniel Killenberger
- 2016 Yusra Alkhazraji
- 2016 Martin Wehrle
- 2014-2015 Patrick von Reth
- 2009-2014 Erez Karpas
- 2014 Robert P. Goldman
- 2010-2012 Andrew Coles
- 2010, 2012 Patrik Haslum
- 2003-2011 Silvia Richter
- 2009-2011 Emil Keyder
- 2010-2011 Moritz Gronbach
- 2010-2011 Manuela Ortlieb
- 2011 Vidal Alcázar Saiz
- 2011 Michael Katz
- 2011 Raz Nissim
- 2010 Moritz Goebelbecker
- 2007-2009 Matthias Westphal
- 2009 Christian Muise


## History

The current version of Fast Downward is the merger of three different
projects:

- the original version of Fast Downward developed by Malte Helmert
  and Silvia Richter
- LAMA, developed by Silvia Richter and Matthias Westphal based on
  the original Fast Downward
- FD-Tech, a modified version of Fast Downward developed by Erez
  Karpas and Michael Katz based on the original code

In addition to these three main sources, the codebase incorporates
code and features from numerous branches of the Fast Downward codebase
developed for various research papers. The main contributors to these
branches are Malte Helmert, Gabi Röger and Silvia Richter.


## License

```
Fast Downward is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

Fast Downward is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```
