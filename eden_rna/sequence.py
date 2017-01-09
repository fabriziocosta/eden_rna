#!/usr/bin/env python
"""Provides trivial sequence folding for RNA sequence."""

from eden_rna import sequence_dotbracket_to_graph


def seq_to_graph(header, sequence):
    """Fold a sequence in a path graph."""
    seq_struct = '.' * len(sequence)
    graph = sequence_dotbracket_to_graph(seq_info=sequence,
                                         seq_struct=seq_struct)
    graph.graph['info'] = 'sequence'
    graph.graph['sequence'] = sequence
    graph.graph['structure'] = seq_struct
    graph.graph['id'] = header
    return graph


def fold(seqs):
    """Fold a list of sequences into path graphs."""
    for header, seq in seqs:
        yield seq_to_graph(header, seq)
