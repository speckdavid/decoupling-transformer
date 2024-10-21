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

import filters

import decoupling_parser

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
BASE_REVISION = "234e359bd7f676af5e40e5ebe3cdd26eb044abce"
REVISION = "cdb75a84f62ac01fc20d4eee3263f1c536add008"
REVISIONS = [BASE_REVISION, REVISION]

CONFIGS = []
factorings = {
    'LP-F0.2s1M-30':        'lp(factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-L0.2s1M-30':        'lp(factoring_time_limit=30, strategy=mml, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-AS0.2s1M-30':       'lp(factoring_time_limit=30, strategy=mmas, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-M0.2s1M-30':        'lp(factoring_time_limit=30, strategy=mm, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-Mo0.2s1M-30':       'lp(factoring_time_limit=30, strategy=mm_opt, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
}
heuristics = {"ff" : [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "[hff], preferred=[hff]"],

}

DRIVER_OPTS = ["--overall-time-limit", "5m"]

for h_name, heuristic_option in heuristics.items():
    h_predefine, search_option = heuristic_option
    CONFIGS.append(IssueConfig(f'{h_name}', h_predefine + ['--search',  f'lazy_greedy({search_option}, cost_type=one)'], driver_options=DRIVER_OPTS))
    for dec_name, dec in factorings.items():
        CONFIGS.append(IssueConfig(f'{h_name}-{dec_name}', h_predefine + ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy({search_option}, cost_type=one)'], driver_options=DRIVER_OPTS))


SUITE = common_setup.DEFAULT_SATISFICING_SUITE

ENVIRONMENT = TetralithEnvironment(
    email="daniel.gnad@liu.se",
#    time_limit_per_task="24:00:00",
#    memory_per_cpu="8300M",
#    extra_options="#SBATCH -A naiss2023-5-236", # parground
    extra_options="#SBATCH -A naiss2023-5-314", # dfsplan
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

def replace_revisions(run):
    run["algorithm"] = run["algorithm"].replace(BASE_REVISION, "base")
    run["algorithm"] = run["algorithm"].replace(REVISION, "fixed")
    return run

exp.add_fetcher(name='fetch', filter=[replace_revisions])


FORMAT = "html"

# REPORT TABLES
attributes = common_setup.ATTRIBUTES

exp.add_report(AbsoluteReport(attributes=attributes), outfile=f"{SCRIPT_NAME}-all.html")

algorithm_pairs = [(f"base-ff-{x}", f"fixed-ff-{x}") for x in factorings.keys()]
exp.add_report(ComparativeReport(attributes=attributes, algorithm_pairs=algorithm_pairs), outfile=f"{SCRIPT_NAME}-compare.html")


# SCATTER PLOTS
def concl_leaf_ratio_as_category(run1, run2):
    if "number_conclusive_leaves" not in run2 or "number_leaf_factors" not in run2:
        return -1
    ratio = 100.0 * run2["number_conclusive_leaves"] / run2["number_leaf_factors"]
    ratio = int(ratio / 20) * 20
    return ratio


PLOT_FORMAT = "png"

for alg1, alg2 in algorithm_pairs:
    for attr in ["transformation_time", "task_size"]:
        exp.add_report(
            ScatterPlotReport(
                attributes=[attr],
                filter_algorithm=[alg1, alg2],
                get_category=lambda r1, r2: r1["domain"],
                format=PLOT_FORMAT,
                show_missing=False,
                relative=True,
            ),
            name=f"plot-{attr.replace('_', '-')}-{alg1}-vs-{alg2}",
        )


exp.run_steps()


