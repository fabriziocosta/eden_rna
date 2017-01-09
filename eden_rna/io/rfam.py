#!/usr/bin/env python
"""Provides io for RNAs."""

import os
import requests
from eden_rna.io.fasta import load as load_fasta


def _save(text, full_out_file_name):
    with open(full_out_file_name, 'w') as f:
        for line in text:
            f.write("%s\n" % line.encode('utf8').strip())


def _rfam_uri(rfam_id):
    # retrieve the seed sequences from RFAM and save them to file
    rfam = 'http://rfam.xfam.org/family/'
    body = '%s/alignment?acc=%s' % (rfam_id, rfam_id)
    fmt = '&format=fastau&download=0'
    uri = rfam + body + fmt
    rfam_dir = 'RNA'
    fname = rfam_id + '.fa'
    if not os.path.exists(rfam_dir):
        os.mkdir(rfam_dir)
    full_out_file_name = os.path.join(rfam_dir, fname)
    if not os.path.isfile(full_out_file_name):
        text = requests.get(uri).text.split('\n')
        _save(text, full_out_file_name)
    return full_out_file_name


def get_rfam_sequence(rfam_id, seq_id=0):
    """Connect to RFAM database and retrieve a single sequence."""
    uri = _rfam_uri(rfam_id)
    seqs = load_fasta(uri)
    for id, (header, seq) in enumerate(seqs):
        if id == seq_id:
            return header, seq


def get_rfam_sequences(rfam_id, seq_ids=None):
    """Connect to RFAM database and retrieve sequences."""
    uri = _rfam_uri(rfam_id)
    seqs = load_fasta(uri)
    if seq_ids is None:
        for header, seq in seqs:
            yield header, seq
    else:
        for id, (header, seq) in enumerate(seqs):
            if id in seq_ids:
                yield header, seq


def load(rfam_id, seq_ids=None):
    """Connect to RFAM database and retrieve sequences."""
    return get_rfam_sequences(rfam_id, seq_ids=seq_ids)
