#!/usr/bin/env python3
"""Download IEDB curated dengue epitopes (T-cell I/II + B-cell) to TSV.

Endpoints (PostgREST):
  https://query-api.iedb.org/tcell_search   per-assay T-cell rows
  https://query-api.iedb.org/bcell_search   per-assay B-cell rows

Filter: parent_source_antigen_source_org_name ILIKE *Dengue*
Pagination: HTTP Range header (PostgREST), max page = 10000.
"""
import csv, json, sys, time, urllib.parse, urllib.request

OUTDIR = '/data/james/ica-dengue/iedb_curated_epitopes'
PAGE = 5000
TIMEOUT = 180
BASE = 'https://query-api.iedb.org'

HDR_OUT = ['epitope_sequence', 'source_protein', 'host_HLA',
           'paper_PMID', 'n_subjects', 'assay_type']

def build_url(endpoint, extra_params):
    q = [('parent_source_antigen_source_org_name', 'ilike.*Dengue*')]
    q.extend(extra_params)
    return f'{BASE}/{endpoint}?' + urllib.parse.urlencode(q)

def fetch(endpoint, extra_params, page=PAGE):
    url = build_url(endpoint, extra_params)
    offset = 0
    while True:
        req = urllib.request.Request(
            url,
            headers={
                'Accept': 'application/json',
                'Range-Unit': 'items',
                'Range': f'{offset}-{offset + page - 1}',
            },
        )
        last_err = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                    body = r.read()
                    data = json.loads(body)
                    cr = r.headers.get('Content-Range', '')
                break
            except Exception as e:
                last_err = e
                if attempt == 3:
                    raise
                time.sleep(5 * (attempt + 1))
                sys.stderr.write(f'retry {endpoint} offset={offset}: {e}\n')
        if not data:
            return
        for row in data:
            yield row
        sys.stderr.write(f'  {endpoint}: {offset + len(data)} rows fetched ({cr})\n')
        if len(data) < page:
            return
        offset += page

def host_hla(row):
    parts = []
    if row.get('mhc_allele_name'):
        parts.append(str(row['mhc_allele_name']))
    if row.get('mhc_restriction') and row.get('mhc_restriction') not in parts:
        parts.append(str(row['mhc_restriction']))
    return ';'.join([p for p in parts if p]) or ''

def source_protein(row):
    csa = row.get('curated_source_antigen') or {}
    if isinstance(csa, dict):
        nm = csa.get('name')
        acc = csa.get('accession')
        if nm and acc:
            return f'{nm} ({acc})'
        if nm:
            return str(nm)
    return str(row.get('parent_source_antigen_name') or '')

def n_subjects(row):
    qm = row.get('qualitative_measure') or ''
    ev = row.get('direct_ex_vivo_bool')
    if ev in (True, 1):
        return f'{qm}; ex_vivo'.strip('; ')
    return qm

def assay_type(row):
    a = row.get('assay_names')
    if isinstance(a, list):
        return '|'.join([str(x) for x in a])
    if a:
        return str(a)
    desc = row.get('assay_description') or ''
    return str(desc).replace('\t', ' ').replace('\n', ' ').replace('<br/>', ' | ')[:300]

def write_tsv(path, endpoint, extra_params):
    n = 0
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        w.writerow(HDR_OUT)
        for row in fetch(endpoint, extra_params):
            seq = row.get('linear_sequence') or row.get('structure_description') or ''
            w.writerow([
                seq,
                source_protein(row),
                host_hla(row),
                row.get('pubmed_id') or '',
                n_subjects(row),
                assay_type(row),
            ])
            n += 1
    return n

def main():
    out = {
        'tcell_class_i.tsv':  ('tcell_search', [('mhc_class', 'eq.I')]),
        'tcell_class_ii.tsv': ('tcell_search', [('mhc_class', 'eq.II')]),
        'bcell.tsv':          ('bcell_search', []),
    }
    summary = {}
    for fname, (ep, params) in out.items():
        path = f'{OUTDIR}/{fname}'
        sys.stderr.write(f'== {fname}: {ep} {params} ==\n')
        n = write_tsv(path, ep, params)
        summary[fname] = n
        sys.stderr.write(f'   -> {n} rows -> {path}\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
