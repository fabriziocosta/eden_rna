#!/usr/bin/env python
"""Provides sequence folding for RNA sequence."""

import os
import subprocess as sp
from eden_rna.sequence import fold as sequence_fold
from eden.util import is_iterable
import math
import networkx as nx
import logging

logger = logging.getLogger(__name__)


def rnaplfold_wrapper(sequence,
                      max_num_edges=None,
                      window_size=None,
                      max_bp_span=None,
                      avg_bp_prob_cutoff=None,
                      no_lonely_bps=None):
    """Fold sequence using RNAplfold."""
    no_lonely_bps_str = ""
    if no_lonely_bps:
        no_lonely_bps_str = "--noLP"
    # Call RNAplfold on command line.
    cmd = 'echo "%s" | RNAplfold -W %d -L %d -c %.2f %s' % (sequence,
                                                            window_size,
                                                            max_bp_span,
                                                            avg_bp_prob_cutoff,
                                                            no_lonely_bps_str)
    sp.check_output(cmd, shell=True)
    # Extract base pair information.
    start_flag = False
    plfold_bp_list = []
    with open('plfold_dp.ps') as f:
        for line in f:
            if start_flag:
                values = line.split()
                if len(values) == 4:
                    avg_prob = values[2]
                    source_id = values[0]
                    dest_id = values[1]
                    plfold_bp_list.append((avg_prob, source_id, dest_id))
            if 'start of base pair probability data' in line:
                start_flag = True
    # Delete RNAplfold output file.
    os.remove("plfold_dp.ps")
    # Return list with base pair information.
    return plfold_bp_list


def _string_to_networkx(header, sequence, **options):
    # Sort edges by average base pair probability in order to stop after
    # max_num_edges edges have been added to a specific vertex.
    max_num_edges = options.get('max_num_edges', 1)
    window_size = options.get('window_size', 150)
    window_size = min(window_size, len(sequence))
    max_bp_span = options.get('max_bp_span', 100)
    max_bp_span = min(max_bp_span, len(sequence))
    avg_bp_prob_cutoff = options.get('avg_bp_prob_cutoff', 0.2)
    no_lonely_bps = options.get('no_lonely_bps', True)
    nesting = options.get('nesting', False)
    hard_threshold = options.get('hard_threshold', False)

    plfold_bp_list = sorted(rnaplfold_wrapper(
        sequence,
        max_num_edges=max_num_edges,
        window_size=window_size,
        max_bp_span=max_bp_span,
        avg_bp_prob_cutoff=avg_bp_prob_cutoff,
        no_lonely_bps=no_lonely_bps), reverse=True)
    graph = nx.Graph()
    graph.graph['id'] = header
    graph.graph['info'] = \
        'RNAplfold: ne=%s ws=%s max_bp_span=%s p_cutoff=%s nolbp=%s' % (
        max_num_edges, window_size, max_bp_span,
        avg_bp_prob_cutoff, no_lonely_bps)
    graph.graph['sequence'] = sequence
    # Add nucleotide vertices.
    for i, c in enumerate(sequence):
        graph.add_node(i, label=c, position=i)
    # Add plfold base pairs and average probabilites.
    for avg_prob_str, source_str, dest_str in plfold_bp_list:
        source = int(source_str) - 1
        dest = int(dest_str) - 1
        avg_prob = math.pow(float(avg_prob_str), 2)
        # Check if either source or dest already have more than
        # max_num_edges edges.
        if len(graph.edges(source)) >= max_num_edges or \
                len(graph.edges(dest)) >= max_num_edges:
            pass
        else:
            if nesting:
                if avg_prob >= hard_threshold:
                    graph.add_edge(source,
                                   dest,
                                   label='=',
                                   type='basepair',
                                   weight=avg_prob)
                else:
                    graph.add_edge(source,
                                   dest,
                                   label='=',
                                   type='basepair',
                                   nesting=True,
                                   weight=avg_prob)
            else:
                graph.add_edge(source,
                               dest,
                               label='=',
                               type='basepair',
                               weight=avg_prob)
    # Add backbone edges.
    for i, c in enumerate(sequence):
        if i > 0:
            graph.add_edge(i, i - 1, label='-', type='backbone')
    return graph


def rnaplfold_to_eden(iterable, **options):
    """Fold RNA sequence with RNAfold."""
    assert(is_iterable(iterable)), 'Not iterable'
    for header, seq in iterable:
        try:
            graph = _string_to_networkx(header, seq, **options)
        except Exception as e:
            logger.debug(e.__doc__)
            logger.debug(e.message)
            logger.debug('Error in: %s' % seq)
            graph = sequence_fold(header, seq, **options)
        yield graph


def fold(seqs, **options):
    """Fold RNA sequence with RNAfold."""
    return rnaplfold_to_eden(seqs, **options)
