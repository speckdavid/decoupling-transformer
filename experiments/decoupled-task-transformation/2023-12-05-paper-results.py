#! /usr/bin/env python3

import itertools
import math
import os
from pathlib import Path
import subprocess

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



attributes = common_setup.ATTRIBUTES

dec_filter = common_setup.NonDecoupledTaskFilter()

exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}.html")

factoring_references = [f"{REVISION}-ff-{name}" for name in [#"L1.0s1M", "L1.0s1M1L", 
                                                             "C1.0s1M", "C1.0s1M1L", #"C0.2s1M", "C0.2s1M1L", 
                                                             "F0.2s1M", "F0.2s1M1L"]]

for ref in factoring_references:
    dec_filter = common_setup.NonDecoupledTaskFilter([ref])

    exp.add_report(AbsoluteReport(attributes=attributes, filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs]), outfile=f"{SCRIPT_NAME}-filter-{ref[41:]}.html")

def filter_only_lama(run):
    if "lama" in run["algorithm"]:
        return run
    return False

def filter_no_lama(run):
    if "lama" not in run["algorithm"]:
        return run
    return False

dec_filter = common_setup.NonDecoupledTaskFilter([f"{REVISION}-ff-F0.2s1M"])
exp.add_report(AbsoluteReport(attributes=["coverage"], filter=[dec_filter.add_runs, dec_filter.filter_non_decoupled_runs, filter_only_lama]), outfile=f"{SCRIPT_NAME}-filter-{ref[41:]}-coverage-lama.html")
exp.add_report(AbsoluteReport(attributes=["coverage"], filter=[dec_filter.filter_non_decoupled_runs, filter_no_lama]), outfile=f"{SCRIPT_NAME}-filter-{ref[41:]}-coverage-ff-po.html")

 

exp.run_steps()

dec_filter.print_statistics()

