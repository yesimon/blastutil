import argparse
import collections
import logging
import os
from os.path import join
import gzip
import time

log = logging.getLogger()

class BlastRecord:
    __slots__ = ('query_id', 'subject_id', 'percent_identity', 'aln_length', 'mismatch_count', 'gap_open_count', 'query_start',
                 'query_end', 'subject_start', 'subject_end', 'e_val', 'bit_score', 'gi', 'accession', 'taxids', 'scinames', 'comnames', 'title',
                 'blast_type', 'line')

    def __init__(self, *args, line=None):
        self.query_id = None
        self.subject_id = None
        self.percent_identity = None
        self.aln_length = None
        self.mismatch_count = None
        self.gap_open_count = None
        self.query_start = None
        self.query_end = None
        self.subject_start = None
        self.subject_end = None
        self.e_val = None
        self.bit_score = None
        self.gi = None
        self.accession = None
        self.taxids = []
        self.scinames = []
        self.comnames = []
        self.title = []
        self.blast_type = None
        for attr, val in zip(self.__slots__, args):
            setattr(self, attr, val)
        self.line = line

    def __eq__(self, other):
        return (self.query_id == other.query_id and self.subject_id == other.subject_id and
                self.query_start == other.query_start and self.query_end == other.query_end and
                self.subject_start == other.subject_start and self.subject_end == other.subject_end)

    def __ne__(self, other):
        """Overrides the default implementation (unnecessary in Python 3)"""
        return not self.__eq__(other)

    def __str__(self):
        if self.line is not None:
            return self.line
        return '\t'.join(str(x) for x in [self.query_id, self.subject_id, self.percent_identity, self.aln_length,
                         self.mismatch_count, self.gap_open_count, self.query_start, self.query_end,
                         self.subject_start, self.subject_end, self.e_val, self.bit_score, self.gi,
                                         self.accession, ';'.join(str(y) for y in self.taxids) if self.taxids else 'N/A', self.scinames,
                                         self.comnames, self.title])

def blast_records(f, blast_type=None):
    '''Yield blast m8 records line by line'''
    for line in f:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        for field in range(3, 10):
            parts[field] = int(parts[field])
        if len(parts) > 12:
            parts[12] = int(parts[12])
            parts[14] = [int(x) for x in parts[14].split(';')
                         if x != 'N/A']
        for field in (2, 10, 11):
            parts[field] = float(parts[field])

        rec = BlastRecord(*parts, line=line)
        if blast_type:
            rec.blast_type = blast_type
        yield rec


def blast_tuples(f, blast_type=None):
    '''Yield blast m8 records line by line'''
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.rstrip().split('\t')
        for field in range(3, 10):
            parts[field] = int(parts[field])
        if len(parts) > 12:
            parts[12] = int(parts[12])
            parts[14] = [int(x) for x in parts[14].split(';')
                         if x != 'N/A']
        for field in (2, 10, 11):
            parts[field] = float(parts[field])

        if blast_type:
            parts.append(blast_type)
        yield parts
