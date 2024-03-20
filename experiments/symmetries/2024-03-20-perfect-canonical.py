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
REVISION = "5e60e606817f953e5ea038cc083dca1e2aa133be"
REVISIONS = [REVISION]

CONFIGS = []
symmetries = ["none", "init", "random"]

searches = {"hmax" :    [[], "astar(hmax())"],
            "blind" : [[], "astar(blind(cache_estimates=false))"]
}

DRIVER_OPTS = ["--overall-memory-limit", "3G", "--overall-time-limit", "5m"]

for h_name, heuristic_options in searches.items():
    h_predefine, search_option = heuristic_options
    for sym_type in symmetries:
        CONFIGS.append(IssueConfig(f'{h_name}-{sym_type[:3]}', h_predefine + ['--root-task-transform', f"symmetry(symmetries=structural_symmetries(), compute_perfect_canonical=false, empty_value_strategy={sym_type})", '--search',  search_option], driver_options=DRIVER_OPTS))
        CONFIGS.append(IssueConfig(f'{h_name}-{sym_type[:3]}PC', h_predefine + ['--root-task-transform', f"symmetry(symmetries=structural_symmetries(), compute_perfect_canonical=true, empty_value_strategy={sym_type})", '--search',  search_option], driver_options=DRIVER_OPTS))


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

exp.add_fetcher("../../../fd-symmetries-full-orbit/experiments/full-orbits/data/2024-03-20-blind-hmax-full-orbit-eval/", name='fetch-oss-baselines', merge=True)
exp.add_fetcher("data/2024-02-19-baselines-eval/", name='fetch-baseline-blind', filter_algorithm=["842fb8778ba9cfc2cce1b886b1670a552e20a548-blind"], merge=True)
exp.add_fetcher("data/2024-03-20-hmax-baseline-eval/", name='fetch-baseline-hmax', merge=True)

attributes = common_setup.ATTRIBUTES + [Attribute("no_symmetries", absolute=True)]

def filter_revision(run):
    run["algorithm"] = run["algorithm"][41:]
    return run

def filter_only_blind(run):
    if "blind" in run["algorithm"]:
        return run
    return False

def filter_only_hmax(run):
    if "hmax" in run["algorithm"]:
        return run
    return False

sym_filter = common_setup.NoSymmetryTaskFilter()
exp.add_report(AbsoluteReport(attributes=attributes, filter=[filter_revision, sym_filter.add_runs, sym_filter.filter_non_symmetry_runs]), outfile=f"{SCRIPT_NAME}.html")
exp.add_report(AbsoluteReport(attributes=attributes, filter=[filter_only_blind, filter_revision, sym_filter.add_runs, sym_filter.filter_non_symmetry_runs]), outfile=f"{SCRIPT_NAME}-blind.html")
exp.add_report(AbsoluteReport(attributes=attributes, filter=[filter_only_hmax, filter_revision, sym_filter.add_runs, sym_filter.filter_non_symmetry_runs]), outfile=f"{SCRIPT_NAME}-hmax.html")


# plots
PLOT_FORMAT = "png"
BASE_REVISION = "87f44d7fa5685421f6b12c6978981491f2e72d1d"
exp.add_report(
    ScatterPlotReport(
        attributes=["expansions_until_last_jump"],
        filter_algorithm=[f"{BASE_REVISION}-blind-oss-hc", f"{BASE_REVISION}-blind-oss-full"],
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-expansions-hc-vs-full",
)


exp.run_steps()

sym_filter.print_statistics()
