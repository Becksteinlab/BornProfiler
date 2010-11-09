#!/usr/bin/env python
"""Setup Born profile calculation with a membrane."""

from membrane import APBSmem

class BP(object):
    def __init__(self, *args, **kwargs):
        self.pqrName = args[0]
        self.pointsName = args[1]
        

        L = APBSmem(self.pqrName, 'L', **kwargs)
