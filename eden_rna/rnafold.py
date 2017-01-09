#!/usr/bin/env python
"""Provides sequence folding for RNA sequence."""


import subprocess as sp
from eden_rna.sequence import fold as sequence_fold
from eden_rna import sequence_dotbracket_to_graph
from eden.util import is_iterable
import logging

logger = logging.getLogger(__name__)


def rnafold_wrapper(sequence, **options):
    """Wrap RNAfold."""
    # defaults
    flags = options.get('flags', '--noPS')
    # command line
    cmd = 'echo "%s" | RNAfold %s' % (sequence, flags)
    out = sp.check_output(cmd, shell=True)
    text = out.strip().split('\n')
    seq_info = text[0]
    seq_struct = text[1].split()[0]
    return seq_info, seq_struct


def _string_to_networkx(header, sequence, **options):
    seq_info, seq_struct = rnafold_wrapper(sequence, **options)
    graph = sequence_dotbracket_to_graph(seq_info=seq_info,
                                         seq_struct=seq_struct)
    graph.graph['info'] = 'RNAfold'
    graph.graph['sequence'] = sequence
    graph.graph['structure'] = seq_struct
    graph.graph['id'] = header
    return graph


def rnafold_to_eden(iterable=None, **options):
    """Fold RNA seq with RNAfold.

    Parameters
    ----------
    iterable: over (header_string, sequence_string)

    options

    Returns
    -------
        nx.graph generator
    """
    assert (is_iterable(iterable)), 'Not iterable'
    for header, seq in iterable:
        try:
            graph = _string_to_networkx(header, seq, **options)
        except Exception as e:
            logger.debug(e.__doc__)
            logger.debug(e.message)
            logger.debug('Error in: %s' % seq)
            graph = sequence_fold(header, seq, **options)
        yield graph


def fold(seqs):
    """Fold RNA sequence with RNAfold."""
    return rnafold_to_eden(seqs)
