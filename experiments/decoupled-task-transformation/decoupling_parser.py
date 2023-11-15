#! /usr/bin/env python

import re

from lab.parser import Parser


class DecouplingParser(Parser):
    def __init__(self):
        super().__init__()
        self.add_pattern('number_leaf_factors', 'Number leaf factors: (.+)', required=False, type=int)
        self.add_pattern('factoring_time', 'Factoring time: (.+)s', required=False, type=float)

        self.add_pattern('ff_simplify_time', 'time to simplify: (.+)s', required=False, type=float)
        self.add_pattern('ff_num_unary_operators', '\[(.+) unary operators\]', required=False, type=int)
        
        self.add_pattern('transformation_time', 'Time for decoupled transformation: (.+)s', required=False, type=float)
        self.add_pattern('number_variables', 'Number of variables: (.+)', required=False, type=int)
        self.add_pattern('number_prime_variables', 'Number of primary variables: (.+)', required=False, type=int)
        self.add_pattern('number_second_variables', 'Number of secondary variables: (.+)', required=False, type=int)
        self.add_pattern('number_operators', 'Number of operators: (.+)', required=False, type=int)
        self.add_pattern('number_axioms', 'Number of axioms: (.+)', required=False, type=int)

