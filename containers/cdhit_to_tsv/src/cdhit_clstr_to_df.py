#!/usr/bin/env python3
"""Convert a CD-HIT .clstr file to a flat TSV table.

PATCH (dengue fork): rewritten in stdlib only. The original used pandas
which triggered an np.bool AttributeError at import inside Pete's
amancevice/pandas:1.4.3 container when apptainer's host-home auto-mount
overlaid the container's numpy. Stdlib avoids the entire issue.
"""
import csv
import sys

cdhit_inpath = sys.argv[1]
cdhit_outpath = sys.argv[2]


def cluster_to_rows(cluster):
    lines = cluster.split('\n')
    index = lines[0].strip(' ')
    rows = []
    for entry in lines[1:]:
        if not entry.strip():
            continue
        entry_items = entry.split(' ')
        entry_index = entry_items[0].split('\t')[0]
        entity = entry_items[1].strip('...')
        if entry_items[-1] == '*':
            is_centroid = "True"
            similarity = "100.0"
        else:
            is_centroid = "False"
            similarity = entry_items[-1].strip('%')
        rows.append([index, entry_index, entity, is_centroid, similarity])
    return rows


with open(cdhit_inpath, 'r') as f:
    raw_data = f.read()

raw_clusters = [c.strip(' ').strip('\n') for c in raw_data.split('>Cluster')[1:]]
raw_clusters = [c for c in raw_clusters if c]

with open(cdhit_outpath, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['cluster_index', 'entry_index', 'entity', 'is_centroid', 'similarity'])
    for cluster in raw_clusters:
        for row in cluster_to_rows(cluster):
            writer.writerow(row)
