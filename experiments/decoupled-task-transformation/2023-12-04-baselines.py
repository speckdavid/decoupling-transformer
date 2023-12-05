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


DRIVER_OPTS = ["--overall-memory-limit", "8G", "--overall-time-limit", "30m"]

CONFIGS = []
CONFIGS.append(IssueConfig("lama-first", [], driver_options=DRIVER_OPTS + ["--alias", 'lama-first']))
CONFIGS.append(IssueConfig('ff-po', ["--evaluator", "hff=ff(transform=adapt_costs(one))", '--search',  f'lazy_greedy([hff], preferred=[hff])'], driver_options=DRIVER_OPTS))
#CONFIGS.append(IssueConfig('gc', ['--search',  f'lazy_greedy([goalcount(cache_estimates=false)])'], driver_options=DRIVER_OPTS))

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

attributes = common_setup.ATTRIBUTES
attributes.extend(exp.DEFAULT_TABLE_ATTRIBUTES)

exp.add_report(AbsoluteReport(attributes=attributes), outfile=f"{SCRIPT_NAME}.html") 

exp.run_steps()

