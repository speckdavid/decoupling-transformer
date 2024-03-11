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
REVISION = "33958478b728695a2d2cd93cdf5c7adfc294bb64"
REVISIONS = [REVISION]

CONFIGS = []
factorings = {
    'F0.2s1M':     'lp(factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'F0.2s1Mnopt': 'lp(factoring_time_limit=30, strategy=mfa, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000), conclusive_leaf_encoding=basic, skip_unnecessary_leaf_effects=false',
    'MF': 'mf()',
}
heuristics = {"ff" : [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "[hff], preferred=[hff]"],
}

DRIVER_OPTS = ["--overall-memory-limit", "8G", "--overall-time-limit", "30m"]

for h_name, heuristic_options in heuristics.items():
    h_predefine, search_option = heuristic_options
    CONFIGS.append(IssueConfig(h_name, h_predefine + ['--search',  f'lazy_greedy({search_option})'], driver_options=DRIVER_OPTS))
    for dec_name, dec in factorings.items():
        CONFIGS.append(IssueConfig(f'{h_name}-{dec_name}', h_predefine + ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy({search_option})'], driver_options=DRIVER_OPTS))

lama_variants = [
    'lama-first',
    'decoupled-lama-first',
    'miura-lama-first',
]

for variant in lama_variants:
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

exp.add_fetcher(name='fetch', filter=[filters.remove_revision])


DEC_REVISION = "243a4d6c4a49e61bb4f442d677891d38b83fc9a8"
exp.add_fetcher("/proj/parground/users/x_dangn/decoupled-fd/experiments/decoupling-transformer/data/2024-03-11-baselines-eval", name="fetch-dec-base", merge=True, filter=[filters.remove_revision, filters.rename_ds_base])

MS_REVISION = "8c2bbe375222e0ba6cdb83b4c79ee4ade6e884ad"
exp.add_fetcher("/proj/parground/users/x_dangn/torralba-sievers-ijcai2019-fast-downward/experiments/decoupling-transformer/data/2024-03-11-ff-po-als-eval", name="fetch-ms", merge=True, filter=[filters.remove_revision])

FORMAT = "tex"

# REPORT TABLES
attributes = common_setup.ATTRIBUTES + [Attribute("compact_coverage", absolute=True, min_wins=False)]

dec_filter = filters.NonDecoupledTaskFilter() # remove instances that are not decoupleable
exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-all.html")

cov_filter = filters.CompactCoverageFilter(["ff-F0.2s1Mnopt", "ff-F0.2s1M", "ff", "DS-ff-F0.2s1M", "RS-a-ls", "decoupled-lama-first", "lama-first"])
dec_filter2 = filters.NonDecoupledTaskFilter(["ff-F0.2s1M"]) # remove instances that are not decoupleable
exp.add_report(AbsoluteReport(attributes=attributes, filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M", "ff", "DS-ff-F0.2s1M", "RS-a-ls", "decoupled-lama-first", "lama-first"], filter=[dec_filter2.add_runs, dec_filter2.filter_non_decoupled_runs, cov_filter.add_run, cov_filter.set_compact_coverage], format=FORMAT), outfile=f"{SCRIPT_NAME}-no-mf.{FORMAT}")

dec_filter4 = filters.NonDecoupledTaskFilter(["ff-MF"]) # remove instances that are not decoupleable
exp.add_report(AbsoluteReport(attributes=attributes, filter_algorithm=["ff-MF", "ff", "miura-lama-first", "lama-first"], filter=[dec_filter4.add_runs, dec_filter4.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-mf.html")




# SCATTER PLOTS
def concl_leaf_ratio_as_category(run1, run2):
    if "number_conclusive_leaves" not in run2 or "number_leaf_factors" not in run2:
        return -1
    ratio = 100.0 * run2["number_conclusive_leaves"] / run2["number_leaf_factors"]
    ratio = int(ratio / 20) * 20
    return ratio


PLOT_FORMAT = "tex"

size_setter = filters.PlotTaskSizeSetter() # hacky way to get the task size of the baseline, which is not printed
exp.add_report(
    ScatterPlotReport(
        attributes=["plot_task_size"],
        filter=[size_setter.get_plot_task_size, size_setter.set_plot_task_size],
        filter_algorithm=["ff", "ff-F0.2s1M"],
        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-transformation",
)
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M"],
        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-optimization",
)




# OTHER STATISTICS
dec_filter3 = filters.NonDecoupledTaskFilter() # remove instances that are not decoupleable
trans_time_check = filters.TranformationTimeChecker("ff-F0.2s1M")
exp.add_report( # dummy report, only collecting statistics
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M"],
        filter=[dec_filter3.add_runs, dec_filter3.filter_non_decoupled_runs, trans_time_check.get_time],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="get-transformation-time-statistics",
)

mf_time_check = filters.MFTimeChecker()
exp.add_report( # dummy report, only collecting statistics
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M"],
        filter=[mf_time_check.get_time],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="get-mf-runtime-statistics",
)
 

exp.run_steps()

dec_filter.print_statistics()

trans_time_check.print_histogram()

mf_time_check.print_statistics()



