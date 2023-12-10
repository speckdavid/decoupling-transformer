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

import decoupling_parser

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
REVISION = "d261f84e10bc1fd8ef4f8e01e92635023e332f34"
REVISIONS = [REVISION]

CONFIGS = []
factorings = {
    'Cs1Lnceffo' : 'lp(factoring_time_limit=30, strategy=mcl, add_cg_sccs=true, min_number_leaves=1)',
    'Cs1M1Lnceffo' : 'lp(factoring_time_limit=30, strategy=mcl, add_cg_sccs=true, max_leaf_size=1000000, min_number_leaves=1)',

    'Ls1Lnceffo' : 'lp(factoring_time_limit=30, strategy=mml, add_cg_sccs=true, min_number_leaves=1)',
    'Ls1M1Lnceffo' : 'lp(factoring_time_limit=30, strategy=mml, add_cg_sccs=true, max_leaf_size=1000000, min_number_leaves=1)',

    'F0.2snceffo': 'lp(factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2)',
    'F0.2s1Mnceffo': 'lp(factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',

    #'MFnceffo': 'mf()',
}
heuristics = {"ff" : [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "[hff], preferred=[hff]"],
              #"gc" : [[], "[goalcount(cache_estimates=false)]"]
}

DRIVER_OPTS = ["--overall-memory-limit", "8G", "--overall-time-limit", "30m"]

for h_name, heuristic_options in heuristics.items():
    h_predefine, search_option = heuristic_options
    for dec_name, dec in factorings.items():
        CONFIGS.append(IssueConfig(f'{h_name}-{dec_name}', h_predefine + ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy({search_option})'], driver_options=DRIVER_OPTS))

lama_variants = [
    'dec-ls1l-lama-first', 
    'dec-ls1m1l-lama-first',
    
    'dec-f02s-lama-first',
    'dec-f02s1m-lama-first',
    
    "dec-cs1l-lama-first",
    "dec-cs1m1l-lama-first",

    #'dec-mf-lama-first',
]

for variant in lama_variants:
    CONFIGS.append(IssueConfig(f"{variant[:-6]}-nceffo", [], driver_options=DRIVER_OPTS + ["--alias", variant]))

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

BASE_REVISION = REVISION 
DEC_REVISION = "613322a5185174645e76b9b63e00673e42e8e95d"
exp.add_fetcher("data/2023-11-30-baselines-5m-eval", name="fetch-expl-base", filter_algorithm=[f"{BASE_REVISION}-ff-po", f"{BASE_REVISION}-gc"], merge=True)
exp.add_fetcher("/proj/parground/users/x_dangn/decoupled-fd/experiments/decoupling-transformer/data/2023-11-22-baselines-5m-eval", name="fetch-dec-base", filter_algorithm=[f"{DEC_REVISION}-ff-L1.0s1M", f"{DEC_REVISION}-ff-F0.2s1M"], merge=True)
exp.add_fetcher("/proj/parground/users/x_dangn/decoupled-fd/experiments/decoupling-transformer/data/2023-11-22-gc-5m-eval", name="fetch-dec-gc", merge=True)

attributes = common_setup.ATTRIBUTES
attributes.extend(exp.DEFAULT_TABLE_ATTRIBUTES)

dec_filter = common_setup.NonDecoupledTaskFilter()

exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}.html")

factoring_references = [f"{REVISION}-ff-{name}" for name in factorings.keys()]

for ref in factoring_references:
    dec_filter = common_setup.NonDecoupledTaskFilter([ref])

    exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-filter-{ref[41:]}.html")


 

exp.run_steps()

dec_filter.print_statistics()

