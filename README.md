# Decoupled Search for the Masses: A Novel Task Transformation for Classical Planning

This repository contains the code for the ICAPS 2024 paper ([[pdf]](https://mrlab.ai/papers/speck-gnad-icaps2024.pdf) [[bib]](https://mrlab.ai/papers/speck-gnad-icaps2024.html)). 
It extends [Fast Downward 23.06](https://github.com/aibasel/downward/) by the possibility of running the search on a decoupled representation of the planning problem.

## Build instructions

See [BUILD.md](BUILD.md). To use the decoupled representation, it is necessary to install CPLEX as described in the file.

## Recommended Search Configurations

We recommend the following search configuration which we defined as aliases.

```
./fast-downward.py --alias XXX domain.pddl problem.pddl
```

See [driver/aliases.py](driver/aliases.py) for more aliases.

## Other Search Configuration

You can also create your own search configuration on the decoupled task as follows.

```
./fast-downward.py --root-task-transform "decoupled(factoring=lp(factoring_time_limit=30, strategy=mcl, add_cg_sccs=true, min_number_leaves=2))" --search "XXX"
```
TODO: use a reasonable config above!
Here, the search "XXX" can be selected as in normal Fast Downward, e.g., `"lazy_greedy([ff()])"`. For more information, see the [Fast Downward website](https://www.fast-downward.org).

### Parameters

```
--root-task-transform decoupled(factoring, same_leaf_preconditons_single_variable=true, conclusive_leaf_encoding=multivalued, skip_unnecessary_leaf_effects=true, dump_task=false, write_sas=false, write_pddl=false)
```

- `factoring` (Factoring): method that computes the factoring (see below)
- `same_leaf_preconditons_single_variable` (bool): The same preconditions of leaves have a single secondary variable.
- `conclusive_leaf_encoding` ({basic, binary, multivalued}): Conclusive leaf encoding.
  - `basic`: no special treatment for conclusive leaves. Operators have conditional effects regarding conclusive leaves.
  - `binary`: primary conclusive leaf variables are represented by binary variables. Operators do not have conditional effects regarding a conclusive leaf; instead, they set the primary variable corresponding to the unique reached leaf state to true and all others to false.
  - `multivalued`: primary conclusive leaf variables are represented using the original variables in a factored manner. Operators do not have conditional effects regarding a conclusive leaf; they simply set the primary leaf variables to the corresponding values of the reached leaf state.
- `skip_unnecessary_leaf_effects` (bool): Skip unnecessary leaf effects for operators that have no influence or are conclusive on the leaf.
- `dump_task` (bool): Dumps the task to the console.
- `write_sas` (bool): Writes the decoupled task to dec_output.sas.
- `write_pddl` (bool): Writes the decoupled task to dec_domain.pddl and dec_problem.pddl.

#### LP Factoring

```
lp(verbosity=normal, min_number_leaves=2, max_leaf_size=infinity, factoring_time_limit=infinity, optimize_leaf_unique_lstate=true, prune_fork_leaf_state_spaces=false, strategy=MML, min_mobility=1, min_flexibility=0, min_fact_flexibility=0, add_cg_sccs=false, max_merge_steps=0, merge_dependent=false, merge_overlapping=false)
```
- `verbosity` ({silent, normal, verbose, debug}): Option to specify the verbosity level.
  - `silent`: only the most basic output
  - `normal`: relevant information to monitor progress
  - `verbose`: full output
  - `debug`: like verbose with additional debug output
- `min_number_leaves` (int): maximum number of leaves
- `max_leaf_size` (int): maximum domain size product of variables in a leaf
- `factoring_time_limit` (int): timeout for computing the factoring
- `optimize_leaf_unique_lstate` (bool): leaves for which every global operator induces a unique leaf state are optimized
- `prune_fork_leaf_state_spaces` (bool): run simulation-based pruning in fork leaves to reduce their state space
- `strategy` ({mml, mmas, mm_opt, mm_approx, mfa, mm, mcl, mcm}): TODO
  - `MML`: maximize mobile leaves
  - `MMAS`: maximize mobile action schemas
  - `MM_OPT`: maximize mobility
  - `MM_APPROX`: maximize mobility (approximation)
  - `MFA`: maximize mobile facts
  - `MM`: maximize mobility (sum)
  - `MCL`: maximize number of mobile conclusive leaves
  - `MCM`: maximize conclusive mobility, i.e. number of conclusive actions (sum)
- `min_mobility` (int): TODO
- `min_flexibility` (double): TODO
- `min_fact_flexibility` (double): TODO
- `add_cg_sccs` (bool): TODO
- `max_merge_steps` (int): TODO
- `merge_dependent` (bool): TODO
- `merge_overlapping` (bool): TODO

#### Miura & Fukunaga factoring

```
mf(verbosity=normal, min_number_leaves=2, max_leaf_size=infinity, factoring_time_limit=infinity, optimize_leaf_unique_lstate=true, prune_fork_leaf_state_spaces=false)
```

- `verbosity` ({silent, normal, verbose, debug}): Option to specify the verbosity level.
 - `silent`: only the most basic output
 - `normal`: relevant information to monitor progress
 - `verbose`: full output
 - `debug`: like verbose with additional debug output
- `min_number_leaves` (int): maximum number of leaves
- `max_leaf_size` (int): maximum domain size product of variables in a leaf
- `factoring_time_limit` (int): timeout for computing the factoring
- `optimize_leaf_unique_lstate` (bool): leaves for which every global operator induces a unique leaf state are optimized
- `prune_fork_leaf_state_spaces` (bool): run simulation-based pruning in fork leaves to reduce their state space

## Decoupled task to SAS or PDDL

TODO: Describe pipeline

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
