#! /usr/bin/env python

import re

from lab.parser import Parser



def parse_factoring_not_possible(content, props):
    props["factoring_possible"] = 1
    props["two_leaf_factoring_possible"] = 1
    
    if re.search("No mobile factoring possible.", content, re.M):
        props["factoring_possible"] = 0

    if re.search("No mobile factoring with at least 2 leaves possible.", content, re.M):
        props["two_leaf_factoring_possible"] = 0

def parse_factoring_type(content, props):
    props["is_lp_factoring"] = 0
    if re.search("Using LP factoring with strategy", content, re.M):
        props["is_lp_factoring"] = 1
    props["is_miura_factoring"] = 0
    if re.search("Using Miura & Fukunaga factoring method.", content, re.M):
        props["is_miura_factoring"] = 1

def parse_exhausted_search_space(content, props):
    props["exhausted_search_space"] = 0
    if re.search("Completely explored state space", content, re.M):
        props["exhausted_search_space"] = 1

class DecouplingParser(Parser):
    def __init__(self):
        Parser.__init__(self)
        self.add_pattern('number_leaf_factors', 'Number leaf factors: (.+)', required=False, type=int)
        self.add_pattern('factoring_time', 'Factoring time: (.+)s', required=False, type=float)

        self.add_pattern('ff_simplify_time', 'time to simplify: (.+)s', required=False, type=float)
        self.add_pattern('ff_num_unary_operators', 'unary operators... done! \[(.+) unary operators\]', required=False, type=int)
        
        self.add_pattern('transformation_time', 'Time for decoupled transformation: (.+)s', required=False, type=float)
        self.add_pattern("task_size", "Task size: (.+)", required=False, type=int)
        self.add_pattern("original_task_size", "Original task size: (.+)", required=False, type=int)
        self.add_pattern('number_variables', 'Number of variables: (.+)', required=False, type=int)
        self.add_pattern('number_prime_variables', 'Number of primary variables: (.+)', required=False, type=int)
        self.add_pattern('number_second_variables', 'Number of secondary variables: (.+)', required=False, type=int)
        self.add_pattern('number_operators', 'Number of operators: (.+)', required=False, type=int)
        self.add_pattern('number_axioms', 'Number of axioms: (.+)', required=False, type=int)

        self.add_pattern('number_conclusive_leaves', "Number of conclusive leaves: (.+)", required=False, type=int)
        self.add_pattern('number_normal_leaves', "Number of normal leaves: (.+)", required=False, type=int)

        self.add_pattern('number_pruned_operators', "Number of pruned operators: (.+)", required=False, type=int)
        
        self.add_pattern('number_wmis_leaf_candidates', "Number final leaf candidates: (.+)", required=False, type=int)

        self.add_function(parse_factoring_not_possible)
        self.add_function(parse_factoring_type)
        self.add_function(parse_exhausted_search_space)
          

