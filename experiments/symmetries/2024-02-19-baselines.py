#! /usr/bin/env python3

import itertools
import math
import os
from pathlib import Path
import subprocess

from lab.environments import TetralithEnvironment
from lab.reports import Attribute, geometric_mean, arithmetic_mean

from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

import common_setup
from common_setup import IssueConfig, IssueExperiment

import symmetry_parser

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
REVISION = "842fb8778ba9cfc2cce1b886b1670a552e20a548"
REVISIONS = [REVISION]

CONFIGS = []

searches = {"ff" :    [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "lazy_greedy([hff], preferred=[hff])"],
            "blind" : [[], "astar(blind(cache_estimates=false))"]
}

DRIVER_OPTS = ["--overall-memory-limit", "3G", "--overall-time-limit", "5m"]

for h_name, heuristic_options in searches.items():
    h_predefine, search_option = heuristic_options
    CONFIGS.append(IssueConfig(f'{h_name}', h_predefine + ['--search',  search_option], driver_options=DRIVER_OPTS))


SUITE = common_setup.DEFAULT_OPTIMAL_SUITE

ENVIRONMENT = TetralithEnvironment(
    email="daniel.gnad@liu.se",
#    time_limit_per_task="24:00:00",
#    memory_per_cpu="8300M",
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
exp.add_parser(symmetry_parser.SymmetryParser())

exp.add_step('build', exp.build)
exp.add_step('start', exp.start_runs)
exp.add_step("parse", exp.parse)
exp.add_fetcher(name='fetch')


attributes = common_setup.ATTRIBUTES
attributes.extend(exp.DEFAULT_TABLE_ATTRIBUTES)

exp.add_report(AbsoluteReport(attributes=attributes), outfile=f"{SCRIPT_NAME}.html")

exp.run_steps()

