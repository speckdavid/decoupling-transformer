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
REVISION = "ed6f486b5025902a97b59cf85833a0009ad932ba"
REVISIONS = [REVISION]

CONFIGS = []
factorings = {
    'WMIS-F0.2s1M':      'wmis(min_number_leaves=1, factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-L0.2s1M':      'wmis(min_number_leaves=1, factoring_time_limit=30, strategy=mml, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-AS0.2s1M':     'wmis(min_number_leaves=1, factoring_time_limit=30, strategy=mmas, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-M0.2s1M':      'wmis(min_number_leaves=1, factoring_time_limit=30, strategy=mm, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-Mo0.2s1M':     'wmis(min_number_leaves=1, factoring_time_limit=30, strategy=mm_opt, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',

}

heuristics = {"ff" : [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "[hff], preferred=[hff]"],
}
DRIVER_OPTS = ["--overall-time-limit", "10m"]

for h_name, heuristic_option in heuristics.items():
    h_predefine, search_option = heuristic_option
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

exp.add_fetcher(name='fetch', filter=[filters.remove_revision])


FORMAT = "html"

# REPORT TABLES
attributes = common_setup.ATTRIBUTES + [Attribute("number_decoupled_tasks", min_wins=False)]

exp.add_fetcher("data/2024-07-24-compare-icaps2024-eval/", name='fetch-lp-factoring', filter_algorithm=[f"new-ff-{x.replace('WMIS', 'LP')}" for x in factorings.keys()], merge=True)

exp.add_report(AbsoluteReport(attributes=attributes, filter=[filters.add_number_decoupled_tasks]), outfile=f"{SCRIPT_NAME}-all.html")

algorithm_pairs = [(f"ff-{x}", f"new-ff-{x.replace('WMIS', 'LP')}") for x in factorings.keys()]
exp.add_report(ComparativeReport(attributes=attributes, algorithm_pairs=algorithm_pairs, filter=[filters.add_number_decoupled_tasks]), outfile=f"{SCRIPT_NAME}-compare.html")

PLOT_FORMAT = "png"

size_setter = filters.PlotTaskSizeSetter() # hacky way to get the task size of the baseline, which is not printed
exp.add_report(
    ScatterPlotReport(
        attributes=["plot_task_size"],
        filter=[size_setter.get_plot_task_size, size_setter.set_plot_task_size],
        filter_algorithm=["ff", "ff-F0.2s1M"],
#        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-transformation",
)
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M"],
#        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-optimization",
)


exp.run_steps()


