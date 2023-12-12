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

#C_ONE_L_REVISION = "32ef3a7297057a014a389c1291c3696f00d22ad4"
#exp.add_fetcher("data/2023-12-05-one-leaf-no-restriction-eval", name="fetch-one-leaf-no-res-C", merge=True, filter=[remove_revision])

#L_ONE_L_REVISION = "e696e5eedaaec073409dbcf42be686437e40b4ad"
#exp.add_fetcher("data/2023-12-06-L-one-leaf-no-restriction-eval", name="fetch-one-leaf-no-res-L", merge=True, filter=[remove_revision])

#NOOPT_REVISION = "32ef3a7297057a014a389c1291c3696f00d22ad4"
#exp.add_fetcher("data/2023-12-06-ff-no-concl-opt-eval", name="fetch-no-concl-opt-one-leaf-no-res-F", merge=True, filter=[remove_revision])

NOOPT_REVISION = "d261f84e10bc1fd8ef4f8e01e92635023e332f34"
exp.add_fetcher("data/2023-12-10-neg-ceff-optimization-eval", name="fetch-neg-ceff-opt", merge=True, filter=[remove_revision])

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
#exp.add_report(AbsoluteReport(attributes=["coverage"], filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs, filter_no_C_factoring, filter_no_lama]), outfile=f"{SCRIPT_NAME}-filter-F0.2s1M-coverage-ff-po.html")

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

BIG_TABLE_ALGORITHMS = ["ff-F0.2s1Mnopt", "ff-F0.2s1Mnceffo", "ff-po", "DS-ff-F0.2s1M", "RS-a-ls", "dec-f02s1m-lama-nceffo", "lama-first"]

def remove_non_final_algorithms(run):
    if run["algorithm"] in BIG_TABLE_ALGORITHMS:
        return run
    return False

FORMAT = "tex"

dec_filter1 = common_setup.NonDecoupledTaskFilter(["ff-F0.2s1Mnceffo"])
comp_coverage1 = CompactCoverageFilter()
exp.add_report(AbsoluteReport(attributes=[Attribute("compact_coverage", absolute=True, min_wins=False)], filter=[remove_non_final_algorithms, dec_filter1.add_runs, dec_filter1.filter_non_decoupled_runs, comp_coverage1.add_run, comp_coverage1.set_compact_coverage], filter_algorithm=BIG_TABLE_ALGORITHMS, format=FORMAT), outfile=f"{SCRIPT_NAME}-filter-F0.2s1M-coverage.{FORMAT}")

MF_TABLE_ALGORITHMS = ["ff-Cs1Lnceffo", "ff-Ls1Lnceffo", "ff-po", "ff-MF", "dec-cs1l-lama-nceffo", "dec-ls1l-lama-nceffo", "lama-first", "dec-mf-lama-first"]

def remove_non_1leaf_algorithms(run):
    if run["algorithm"] in MF_TABLE_ALGORITHMS:
        return run
    return False

dec_filter2 = common_setup.NonDecoupledTaskFilter(["ff-MF"])
comp_coverage2 = CompactCoverageFilter()
exp.add_report(AbsoluteReport(attributes=[Attribute("compact_coverage", absolute=True, min_wins=False)], filter=[remove_non_1leaf_algorithms, dec_filter2.add_runs, dec_filter2.filter_non_decoupled_runs, comp_coverage2.add_run, comp_coverage2.set_compact_coverage], filter_algorithm=MF_TABLE_ALGORITHMS, format=FORMAT), outfile=f"{SCRIPT_NAME}-filter-mf-coverage-1L.{FORMAT}")




class PlotTaskSizeSetter:
    def __init__(self):
        self.original_sizes = defaultdict(int)
    
    def set_plot_task_size(self, run):
        size = None
        inst = f"{run['domain']}-{run['problem']}"
        if run["algorithm"] == "ff-po" and inst in self.original_sizes:
            size = self.original_sizes[inst]
        elif run["algorithm"] == "ff-F0.2s1Mnceffo" and "task_size" in run:
            size = run["task_size"]
        run["plot_task_size"] = size
        return run

    def get_plot_task_size(self, run):
        if run["algorithm"] == "ff-F0.2s1Mnceffo" and "original_task_size" in run:
            self.original_sizes[f"{run['domain']}-{run['problem']}"] = run["original_task_size"]
        return run

def concl_leaf_ratio_as_category(run1, run2):
    if "number_conclusive_leaves" not in run2 or "number_leaf_factors" not in run2:
        return -1
    ratio = 100.0 * run2["number_conclusive_leaves"] / run2["number_leaf_factors"]
    ratio = int(ratio / 20) * 20
    return ratio


PLOT_FORMAT = "tex"

size_setter = PlotTaskSizeSetter()
exp.add_report(
    ScatterPlotReport(
        attributes=["plot_task_size"],
        filter=[size_setter.get_plot_task_size, size_setter.set_plot_task_size],
        filter_algorithm=["ff-po", "ff-F0.2s1Mnceffo"],
        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-transformation",
)
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1Mnceffo"],
        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-optimization",
)
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1M", "ff-F0.2s1Mnceffo"],
        get_category=concl_leaf_ratio_as_category,
        format=PLOT_FORMAT,
        show_missing=False,
    ),
    name="scatterplot-task-size-nceffo",
)

class TranformationTimeChecker:
    def __init__(self, config):
        self.config = config
        self.times = [0] * 7
        self.max_time = 0
    def get_time(self, run):
        if run["algorithm"] == self.config:
            if "transformation_time" in run:
                time = run["transformation_time"]
                self.max_time = max(time, self.max_time)
                if time < 1:
                    self.times[0] += 1
                elif time < 5:
                    self.times[1] += 1
                elif time < 10:
                    self.times[2] += 1
                elif time < 30:
                    self.times[3] += 1
                elif time < 60:
                    self.times[4] += 1
                else:
                    self.times[5] += 1
            else:
                if "number_leaf_factors" in run:
                    self.times[6] += 1
                else:
                    print(f"unknown: {run['domain']}:{run['problem']}")
        return run
    def print_histogram(self):
        print(f"Transformation time statistics: max={self.max_time}s")
        print("<1s\t<5s\t<10s\t<30s\t<60s\t>=60s\tDNF")
        print("\t".join(str(x) for x in self.times))

dec_filter3 = common_setup.NonDecoupledTaskFilter(["ff-F0.2s1Mnceffo"])
trans_time_check = TranformationTimeChecker("ff-F0.2s1Mnceffo")
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1Mnceffo"],
        filter=[dec_filter3.add_runs, dec_filter3.filter_non_decoupled_runs, trans_time_check.get_time],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="get-transformation-time-statistics",
)

class MFTimeChecker:
    def __init__(self):
        self.times = defaultdict(lambda: list([-1, -1]))
    def get_time(self, run):
        if run["coverage"] == 1 and run["algorithm"] in ["ff-MF", "ff-po", "dec-mf-lama-first", "lama-first"]:
            run_id = f"{run['domain']}:{run['problem']}"
            if run["algorithm"] == "ff-MF":
                self.times[f"ff:{run_id}"][1] = run["planner_time"]
            elif run["algorithm"] == "ff-po":
                self.times[f"ff:{run_id}"][0] = run["planner_time"]
            elif run["algorithm"] == "dec-mf-lama-first":
                self.times[f"lama:{run_id}"][1] = run["planner_time"]
            elif run["algorithm"] == "lama-first":
                self.times[f"lama:{run_id}"][0] = run["planner_time"]
        return run
    def print_statistics(self):
        if self.times.values():
            entries = list(self.times.values())
            assert all(len(x) == 2 for x in entries)
            ratios = [x[0] / x[1] for x in entries if all(e >= 0 for e in x)]
            print(f"Max speedup of MF over SAS baseline:  {max(ratios)}")
            print(f"geometric mean speedup of MF over SAS baseline:   {geometric_mean(ratios)}")
            print(f"artithmetic mean speedup of MF over SAS baseline: {arithmetic_mean(ratios)}")
            print(len(ratios))

mf_time_check = MFTimeChecker()
exp.add_report(
    ScatterPlotReport(
        attributes=["task_size"],
        filter_algorithm=["ff-F0.2s1Mnopt", "ff-F0.2s1Mnceffo"],
        filter=[mf_time_check.get_time],
        #get_category=domain_as_category,
        format="png",  # Use "tex" for pgfplots output.
    ),
    name="get-mf-runtime-statistics",
)


exp.run_steps()

trans_time_check.print_histogram()

mf_time_check.print_statistics()

