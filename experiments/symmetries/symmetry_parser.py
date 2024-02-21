#! /usr/bin/env python

import re

from lab.parser import Parser


def parse_exhausted_search_space(content, props):
    props["exhausted_search_space"] = 0
    if re.search("Completely explored state space", content, re.M):
        props["exhausted_search_space"] = 1

def parse_no_symmetries(content, props):
    props["no_symmetries"] = 0
    if re.search("No symmetries found, aborting..", content, re.M):
        props["no_symmetries"] = 1

        

class SymmetryParser(Parser):
    def __init__(self):
        Parser.__init__(self)
        
        self.add_pattern('transformation_time', 'Time for symmetry transformation: (.+)s', required=False, type=float)
        self.add_pattern("task_size", "Task size: (.+)", required=False, type=int)
        self.add_pattern("original_task_size", "Original task size: (.+)", required=False, type=int)

        self.add_pattern("number_generators", "Total number of generators: (.+)", required=False, type=int)

        self.add_function(parse_exhausted_search_space)
        self.add_function(parse_no_symmetries)
         

