#! /usr/bin/env python3

import itertools
import math
import os
from pathlib import Path
import subprocess
from collections import defaultdict

from lab.reports import Attribute, geometric_mean, arithmetic_mean

from lab.experiment import Experiment

from downward.reports.absolute import AbsoluteReport
from downward.reports.compare import ComparativeReport
from downward.reports.scatter import ScatterPlotReport

import common_setup

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]

exp = Experiment()

def remove_revision(run):
    run["algorithm"] = run["algorithm"][41:]
    return run

def rename_ds_base(run):
    run["algorithm"] = "DS-" + run["algorithm"]
    return run

REVISION = "2df19052f3f65488fa3d4db24cae67b81d3dab93"
exp.add_fetcher("data/2023-12-05-lama-eval", name="fetch-lama", filter=[remove_revision])
exp.add_fetcher("data/2023-12-05-ff-eval", name="fetch-ff", merge=True, filter=[remove_revision])

exp.add_fetcher("data/2023-12-04-baselines-eval", name="fetch-expl-base", merge=True, filter=[remove_revision])

DEC_REVISION = "243a4d6c4a49e61bb4f442d677891d38b83fc9a8"
exp.add_fetcher("/proj/parground/users/x_dangn/decoupled-fd/experiments/decoupling-transformer/data/2023-12-05-baselines-eval", name="fetch-dec-base", merge=True, filter=[remove_revision, rename_ds_base])

MS_REVISION = "8c2bbe375222e0ba6cdb83b4c79ee4ade6e884ad"
exp.add_fetcher("/proj/parground/users/x_dangn/torralba-sievers-ijcai2019-fast-downward/experiments/decoupling-transformer/data/2023-12-04-ff-po-als-eval", name="fetch-ms", merge=True, filter=[remove_revision])

#ONE_L_REVISION = "32ef3a7297057a014a389c1291c3696f00d22ad4"
#exp.add_fetcher("data/2023-12-05-one-leaf-no-restriction-eval", name="fetch-one-leaf-no-res", merge=True, filter=[remove_revision])

NOOPT_ONE_L_REVISION = "32ef3a7297057a014a389c1291c3696f00d22ad4"
exp.add_fetcher("data/2023-12-06-ff-no-concl-opt-eval", name="fetch-no-concl-opt-single-leaf-F", merge=True, filter=[remove_revision])


attributes = common_setup.ATTRIBUTES

dec_filter = common_setup.NonDecoupledTaskFilter()

exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}.html")

#factoring_references = [f"ff-{name}" for name in [#"L1.0s1M", "L1.0s1M1L", 
#                                                  "C1.0s1M", "C1.0s1M1L", #"C0.2s1M", "C0.2s1M1L", 
#                                                  "F0.2s1M", "F0.2s1M1L"]]
#
#for ref in factoring_references:
#    dec_filter = common_setup.NonDecoupledTaskFilter([ref])
#    exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-filter-{ref[3:]}.html")
#
#def filter_only_lama(run):
#    if "lama" in run["algorithm"]:
#        return run
#    return False
#
#def filter_no_lama(run):
#    if "lama" not in run["algorithm"]:
#        return run
#    return False
#
#def filter_no_C_factoring(run):
#    if "C1.0" in run["algorithm"] or "c1" in run["algorithm"]:
#        return False
#    return run
#
#dec_filter = common_setup.NonDecoupledTaskFilter([f"ff-F0.2s1M"])
#exp.add_report(AbsoluteReport(attributes=["coverage"], filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs, filter_no_C_factoring, filter_only_lama]), outfile=f"{SCRIPT_NAME}-filter-F0.2s1M-coverage-lama.html")
#exp.add_report(AbsoluteReport(attributes=["coverage"], filter=[dec_filter.filter_non_decoupled_runs, filter_no_C_factoring, filter_no_lama]), outfile=f"{SCRIPT_NAME}-filter-F0.2s1M-coverage-ff-po.html")

class CompactCoverageFilter():
    def __init__(self):
        self.coverage = defaultdict(lambda : defaultdict(int))
    def add_run(self, run):
        if run["coverage"] == 1:
            self.coverage[run["algorithm"]][run["domain"]] += 1
        return run
    def set_compact_coverage(self, run):
        run["compact_coverage"] = run["coverage"]
        if all(self.coverage[list(self.coverage.keys())[0]][run["domain"]] == self.coverage[alg][run["domain"]] for alg in self.coverage.keys()):
            run["problem"] = run["domain"] + run["problem"]
            run["domain"] = "zzOther"
        return run

BIG_TABLE_ALGORITHMS = ["ff-F0.2s1Mnopt", "ff-F0.2s1M", "ff-po", "DS-ff-F0.2s1M", "RS-a-ls", "dec-f02s1m-lama-first", "lama-first"]

def remove_non_final_algorithms(run):
    if run["algorithm"] in BIG_TABLE_ALGORITHMS:
        return run
    return False

FORMAT = "tex"

dec_filter = common_setup.NonDecoupledTaskFilter(["ff-F0.2s1M"])
comp_coverage = CompactCoverageFilter()
exp.add_report(AbsoluteReport(attributes=[Attribute("compact_coverage", absolute=True, min_wins=False)], filter=[remove_non_final_algorithms, dec_filter.add_runs, dec_filter.filter_non_decoupled_runs, comp_coverage.add_run, comp_coverage.set_compact_coverage], filter_algorithm=BIG_TABLE_ALGORITHMS, format=FORMAT), outfile=f"{SCRIPT_NAME}-filter-F0.2s1M-coverage.{FORMAT}")

MF_TABLE_ALGORITHMS = ["ff-Fs1L", "ff-po", "ff-MF", "dec-cs1l-lama-first", "lama-first", "dec-mf-lama-first"]

def remove_non_1leaf_algorithms(run):
    if run["algorithm"] in MF_TABLE_ALGORITHMS:
        return run
    return False

dec_filter = common_setup.NonDecoupledTaskFilter(["ff-MF"])
comp_coverage = CompactCoverageFilter()
exp.add_report(AbsoluteReport(attributes=[Attribute("compact_coverage", absolute=True, min_wins=False)], filter=[remove_non_1leaf_algorithms, dec_filter.add_runs, dec_filter.filter_non_decoupled_runs, comp_coverage.add_run, comp_coverage.set_compact_coverage], filter_algorithm=MF_TABLE_ALGORITHMS, format=FORMAT), outfile=f"{SCRIPT_NAME}-filter-mf-coverage-1L.{FORMAT}")




class PlotTaskSizeSetter:
    def __init__(self):
        self.original_sizes = defaultdict(int)
    
    def set_plot_task_size(self, run):
        size = None
        inst = f"{run['domain']}-{run['problem']}"
        if run["algorithm"] == "ff-po" and inst in self.original_sizes:
            size = self.original_sizes[inst]
        elif run["algorithm"] == "ff-F0.2s1M" and "task_size" in run:
            size = run["task_size"]
        run["plot_task_size"] = size
        return run

    def get_plot_task_size(self, run):
        if run["algorithm"] == "ff-F0.2s1M" and "original_task_size" in run:
            self.original_sizes[f"{run['domain']}-{run['problem']}"] = run["original_task_size"]
        return run

size_setter = PlotTaskSizeSetter()
exp.add_report(
    ScatterPlotReport(
        attributes=["plot_task_size"],
        filter=[size_setter.get_plot_task_size, size_setter.set_plot_task_size],
        filter_algorithm=["ff-po", "ff-F0.2s1M"],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="scatterplot-task-size-transformation",
)
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1M"],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="scatterplot-task-size-optimization",
)


exp.run_steps()

dec_filter.print_statistics()

