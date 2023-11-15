#! /usr/bin/env python3

import itertools
import math
import os
from pathlib import Path
import subprocess

from lab.environments import LocalEnvironment, TetralithEnvironment
from lab.reports import Attribute, geometric_mean, arithmetic_mean

from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

import common_setup
from common_setup import IssueConfig, IssueExperiment

import decoupling_parser

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
REVISION = "23cfc58a817599f5766a1cba3408f89f264c534c"
REVISIONS = [REVISION]

CONFIGS = []
factorings = {
    'L1.0s' : 'lp(factoring_time_limit=30, strategy=mml, min_flexibility=1.0, add_cg_sccs=true)',
    'L0.2' : 'lp(factoring_time_limit=30, strategy=mml, min_flexibility=0.2)',
    'L' : 'lp(factoring_time_limit=30, strategy=mml)',
}
heuristics = {"ff" : "ff(transform=adapt_costs(one))",
              "gc" : "goalcount(cache_estimates=false)"}


DRIVER_OPTS = ["--overall-memory-limit", "8G", "--overall-time-limit", "5m"]

for h_name, heuristic in heuristics.items():
    CONFIGS.append(IssueConfig(f'{h_name}', ['--search',  f'lazy_greedy([{heuristic}])'], driver_options=DRIVER_OPTS))
    for dec_name, dec in factorings.items():
        CONFIGS.append(IssueConfig(f'{dec_name}-{h_name}', ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy([{heuristic}])'], driver_options=DRIVER_OPTS))
        CONFIGS.append(IssueConfig(f'{dec_name}-{h_name}-deb', ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy([{heuristic}])'], build_options=['debug'], driver_options=DRIVER_OPTS + ['--build', 'debug']))


SUITE = common_setup.DEFAULT_SATISFICING_SUITE

ENVIRONMENT = TetralithEnvironment(
    email="daniel.gnad@liu.se",
#    time_limit_per_task="24:00:00",
    memory_per_cpu="8300M",
    extra_options="#SBATCH -A naiss2023-5-236", # parground
)

exp = IssueExperiment(
    revisions=REVISIONS,
    configs=CONFIGS,
    environment=ENVIRONMENT,
)

exp.add_suite(BENCHMARKS_DIR, SUITE)

exp.add_parser(exp.EXITCODE_PARSER)
exp.add_parser(exp.TRANSLATOR_PARSER)
exp.add_parser(exp.SINGLE_SEARCH_PARSER)
exp.add_parser(exp.PLANNER_PARSER)
exp.add_parser(decoupling_parser.DecouplingParser())

exp.add_step('build', exp.build)
exp.add_step('start', exp.start_runs)
exp.add_step("parse", exp.parse)
exp.add_fetcher(name='fetch')

attributes=[
    Attribute('search_out_of_memory', absolute=True, min_wins=True),
    Attribute('search_out_of_time', absolute=True, min_wins=True),

    # decoupled attributes
    Attribute('factoring_time', absolute=False, min_wins=True, function=arithmetic_mean),
    Attribute("number_leaf_factors", absolute=False),

    # general search attributes
    Attribute('reopened_until_last_jump', absolute=False, min_wins=True, function=arithmetic_mean),
    Attribute('inconsistent_heuristic', absolute=True, min_wins=True, function=sum),
]

attributes.extend(exp.DEFAULT_TABLE_ATTRIBUTES)

exp.add_absolute_report_step(attributes=attributes)

 

exp.run_steps()
