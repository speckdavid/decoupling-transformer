from collections import defaultdict
from lab.reports import geometric_mean, arithmetic_mean

class NonDecoupledTaskFilter:
    def __init__(self, reference_configs=None):
        self.decoupled_tasks = defaultdict(set)
        self.factoring_not_possible_tasks = defaultdict(set)
        self.maybe_not_possible_tasks = defaultdict(set)

        self.unsupported_tasks = defaultdict(set)
        self.translate_oom_tasks = defaultdict(set)

        self.reference_configs = reference_configs

    def add_runs(self, run):
        domain = run["domain"]
        problem = run["problem"]
        if run["error"] == "search-unsupported":
            self.unsupported_tasks[domain].add(problem)
        if run["error"] == "translate-out-of-memory":
            self.translate_oom_tasks[domain].add(problem)
        if self.reference_configs and run["algorithm"] not in self.reference_configs:
            return run
        if "number_leaf_factors" in run and run["number_leaf_factors"] > 0:
            self.decoupled_tasks[domain].add(problem)
        elif ("is_lp_factoring" in run and run["is_lp_factoring"] == 1) or ("is_miura_factoring" in run and run["is_miura_factoring"] == 1):
            self.maybe_not_possible_tasks[domain].add(problem)
        if "factoring_possible" in run and run["factoring_possible"] == 0:
            self.factoring_not_possible_tasks[domain].add(problem) 
        return run

    def filter_non_decoupled_runs(self, run):
        problem = run["problem"]
        domain = run["domain"]
        if problem in self.unsupported_tasks[domain]:
            return False 
        if problem in self.translate_oom_tasks[domain]:
            return False
        if problem in self.factoring_not_possible_tasks[domain]:
            return False
        if problem in self.maybe_not_possible_tasks[domain] and problem not in self.decoupled_tasks[domain]:
            return False
        return run

    def print_statistics(self):
        print(f"Number instances where factoring is not possible: {sum(len(tasks) for tasks in self.factoring_not_possible_tasks.values())}")
        print(f"Number instances where factoring was found: {sum(len(tasks) for tasks in self.decoupled_tasks.values())}")
        print(f"Number instances where factoring may not be possible: {sum(len(tasks) for tasks in self.maybe_not_possible_tasks.values())}")
        print(f"Number unsupported instances: {sum(len(tasks) for tasks in self.unsupported_tasks.values())}")


class MFTimeChecker:
    def __init__(self):
        self.times = defaultdict(lambda: list([-1, -1]))
    def get_time(self, run):
        if run["coverage"] == 1 and run["algorithm"] in ["ff-MF", "ff", "miura-lama-first", "lama-first"]:
            run_id = f"{run['domain']}:{run['problem']}"
            if run["algorithm"] == "ff-MF":
                self.times[f"ff:{run_id}"][1] = run["planner_time"]
            elif run["algorithm"] == "ff":
                self.times[f"ff:{run_id}"][0] = run["planner_time"]
            elif run["algorithm"] == "miura-lama-first":
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


class PlotTaskSizeSetter:
    def __init__(self):
        self.original_sizes = defaultdict(int)
    
    def set_plot_task_size(self, run):
        size = None
        inst = f"{run['domain']}-{run['problem']}"
        if run["algorithm"] == "ff" and inst in self.original_sizes:
            size = self.original_sizes[inst]
        elif run["algorithm"] == "ff-F0.2s1M" and "task_size" in run:
            size = run["task_size"]
        run["plot_task_size"] = size
        return run

    def get_plot_task_size(self, run):
        if run["algorithm"] == "ff-F0.2s1M" and "original_task_size" in run:
            self.original_sizes[f"{run['domain']}-{run['problem']}"] = run["original_task_size"]
        return run

class CompactCoverageFilter():
    def __init__(self, configs):
        self.coverage = defaultdict(lambda : defaultdict(int))
        self.configs = configs
    def add_run(self, run):
        if run["algorithm"] not in self.configs:
            return run
        if run["coverage"] == 1:
            self.coverage[run["algorithm"]][run["domain"]] += 1
        return run
    def set_compact_coverage(self, run):
        run["compact_coverage"] = run["coverage"]
        if all(self.coverage[list(self.coverage.keys())[0]][run["domain"]] == self.coverage[alg][run["domain"]] for alg in self.coverage.keys()):
            run["problem"] = run["domain"] + run["problem"]
            run["domain"] = "zzOther"
        return run

def remove_revision(run):
    run["algorithm"] = run["algorithm"][41:]
    return run

def rename_ds_base(run):
    run["algorithm"] = "DS-" + run["algorithm"]
    return run

def remove_mf_configs(run):
    if any(x in run["algorithm"] for x in ["MF", "miura"]):
        return False
    return run

def add_number_decoupled_tasks(run):
    run["number_decoupled_tasks"] = 0
    if "number_leaf_factors" in run:
        run["number_decoupled_tasks"] = 1
    return run

