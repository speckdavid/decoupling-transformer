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
REVISION = "c2ebc938e69ce3ad620d7eef9b88704a1fd6888f"
REVISIONS = [REVISION]

CONFIGS = []
symmetries = ["none", "init", "random"]

searches = {#"ff" :    [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "lazy_greedy([hff], preferred=[hff])"],
            "blind" : [[], "astar(blind(cache_estimates=false))"]
}

DRIVER_OPTS = ["--overall-memory-limit", "3G", "--overall-time-limit", "5m"]

for h_name, heuristic_options in searches.items():
    h_predefine, search_option = heuristic_options
    for sym_type in symmetries:
        CONFIGS.append(IssueConfig(f'{h_name}-{sym_type}', h_predefine + ['--root-task-transform', f"symmetry(symmetries=structural_symmetries(), empty_value_strategy={sym_type})", '--search',  search_option], driver_options=DRIVER_OPTS))


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

exp.add_fetcher("../../../fd-symmetries-full-orbit/experiments/full-orbits/data/2024-02-15-blind-full-orbits-eval/", name='fetch-oss-blind', filter_algorithm=["d823945b7e5ba223a9ef5f44a5952abf7c3c20f4-blind-oss-hc"], merge=True)
#exp.add_fetcher("../../../fd-symmetries-full-orbit/experiments/full-orbits/data/2024-02-20-ff-full-orbits-eval/", name='fetch-oss-ff', filter_algorithm=["d823945b7e5ba223a9ef5f44a5952abf7c3c20f4-ff-oss-hc"], merge=True)
exp.add_fetcher("data/2024-02-19-baselines-eval/", name='fetch-baselines', merge=True)

attributes = common_setup.ATTRIBUTES + [Attribute("no_symmetries", absolute=True)]

def filter_revision(run):
    run["algorithm"] = run["algorithm"][41:]
    if "ff" in run["algorithm"]:
        return False
    return run

sym_filter = common_setup.NoSymmetryTaskFilter()

exp.add_report(AbsoluteReport(attributes=attributes, filter=[filter_revision, sym_filter.add_runs, sym_filter.filter_non_symmetry_runs]), outfile=f"{SCRIPT_NAME}.html")

exp.run_steps()

