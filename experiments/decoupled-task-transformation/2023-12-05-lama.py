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
REVISION = "2df19052f3f65488fa3d4db24cae67b81d3dab93"
REVISIONS = [REVISION]

CONFIGS = []
variants = [
#    'dec-l1s1m-lama-first',
#    'dec-l1s1m1l-lama-first',
    
    'dec-f02s1m-lama-first',
    'dec-f02s1m1l-lama-first',
    
    "dec-c1s1m-lama-first",
    "dec-c1s1m1l-lama-first",
#    "dec-c02s1m-lama-first",
#    "dec-c02s1m1l-lama-first",

    'dec-mf-lama-first',
]


DRIVER_OPTS = ["--overall-memory-limit", "8G", "--overall-time-limit", "30m"]

for variant in variants:
    CONFIGS.append(IssueConfig(variant, [], driver_options=DRIVER_OPTS + ["--alias", variant]))


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

BASE_REVISION = "1593b8a55827565c7f7d1259f324c82d242ae9f6"
exp.add_fetcher("data/2023-11-22-baselines-5m-eval", filter_algorithm=[f"{BASE_REVISION}-lama-first"], merge=True)

attributes = common_setup.ATTRIBUTES
attributes.extend(exp.DEFAULT_TABLE_ATTRIBUTES)

dec_filter = common_setup.NonDecoupledTaskFilter()

exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}.html")

factoring_references = [#f'{REVISION}-dec-l1s-lama-first',
                        f'{REVISION}-dec-l1s1m-lama-first',
                        #f'{REVISION}-dec-f02s-lama-first',
                        f'{REVISION}-dec-f02s1m-lama-first',
                        f'{REVISION}-dec-mf-lama-first',]

for ref in factoring_references:
    dec_filter = common_setup.NonDecoupledTaskFilter([ref])

    exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-filter-{ref[41:]}.html")


 

exp.run_steps()

dec_filter.print_statistics()

