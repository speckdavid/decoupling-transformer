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
BASE_REVISION = "icaps2024"
REVISION = "9135602a7e2fdd2aa752434a5fc85ccca6668327"
REVISIONS = [REVISION]

CONFIGS = []
factorings = {
    'WMIS-F0.2s1M':        'wmis(factoring_time_limit=30, strategy=mfa,    add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-L0.2s1M':        'wmis(factoring_time_limit=30, strategy=mml,    add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-AS0.2s1M':       'wmis(factoring_time_limit=30, strategy=mmas,   add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-M0.2s1M':        'wmis(factoring_time_limit=30, strategy=mm,     add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'WMIS-Ma0.2s1M':       'wmis(factoring_time_limit=30, strategy=mm_opt, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',

    'LP-F0.2s1M':        'lp(factoring_time_limit=30, strategy=mfa,       add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-L0.2s1M':        'lp(factoring_time_limit=30, strategy=mml,       add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-AS0.2s1M':       'lp(factoring_time_limit=30, strategy=mmas,      add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-M0.2s1M':        'lp(factoring_time_limit=30, strategy=mm,        add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',
    'LP-Ma0.2s1M':       'lp(factoring_time_limit=30, strategy=mm_approx, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000)',

}
heuristics = {"ff" : [["--evaluator", "hff=ff(transform=adapt_costs(one))"], "[hff], preferred=[hff]"],

}

DRIVER_OPTS = ["--overall-time-limit", "10m"]

for h_name, heuristic_option in heuristics.items():
    h_predefine, search_option = heuristic_option
    for dec_name, dec in factorings.items():
        CONFIGS.append(IssueConfig(f'{h_name}-{dec_name}', h_predefine + ['--root-task-transform', f"decoupled(factoring={dec})", '--search',  f'lazy_greedy({search_option}, cost_type=one)'], driver_options=DRIVER_OPTS))

CONFIGS.append(IssueConfig(f'dec-lama', [], driver_options=DRIVER_OPTS + ["--alias", "decoupled-lama-first"]))
CONFIGS.append(IssueConfig(f'dec-lama-wmis', ['--root-task-transform', "decoupled(factoring=wmis(factoring_time_limit=30, strategy=mml, add_cg_sccs=true, min_flexibility=0.2, max_leaf_size=1000000))", "--search", "let(hlm, landmark_sum(lm_factory=lm_reasonable_orders_hps(lm_rhw()),transform=adapt_costs(one),pref=false),let(hff, ff(transform=adapt_costs(one)),lazy_greedy([hff,hlm],preferred=[hff,hlm],cost_type=one,reopen_closed=false)))"], driver_options=DRIVER_OPTS))

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
    run["algorithm"] = run["algorithm"].replace(REVISION, "new")
    return run

exp.add_fetcher(name='fetch', filter=[replace_revisions])

def only_base_revision(run):
    if (run["algorithm"].startswith("base")):
        return run
    return False

exp.add_fetcher("data/2024-07-24-compare-icaps2024-eval/", name='fetch-base-revision', filter=[only_base_revision], merge=True)


FORMAT = "html"

# REPORT TABLES
attributes = common_setup.ATTRIBUTES

def all_but_base_revision(run):
    if (run["algorithm"].startswith("base")):
        return False
    return run


exp.add_report(AbsoluteReport(attributes=attributes, filter=[all_but_base_revision]), outfile=f"{SCRIPT_NAME}-all.html")

algorithm_pairs = [(f"base-ff-{x}", f"new-ff-{x}") for x in factorings.keys() if "LP" in x] + [("base-dec-lama", "new-dec-lama")]
exp.add_report(ComparativeReport(attributes=attributes, algorithm_pairs=algorithm_pairs, filter=[filters.add_number_decoupled_tasks]), outfile=f"{SCRIPT_NAME}-compare-icaps2024.html")

algorithm_pairs = [(f"new-ff-{x}", f"new-ff-{x.replace('LP', 'WMIS')}") for x in factorings.keys() if "LP" in x] + [("new-dec-lama", "new-dec-lama-wmis")]
exp.add_report(ComparativeReport(attributes=attributes, algorithm_pairs=algorithm_pairs), outfile=f"{SCRIPT_NAME}-compare-wmis.html")



# SCATTER PLOTS
PLOT_FORMAT = "png"

algorithm_pairs = [(f"base-ff-{x}", f"new-ff-{x}") for x in factorings.keys() if "LP" in x] + [("base-dec-lama", "new-dec-lama")]

for alg1, alg2 in algorithm_pairs:
    for attr in ["transformation_time", "task_size", "factoring_time", "number_leaf_factors"]:
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


