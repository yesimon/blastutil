#!/usr/bin/env python3
import argparse
import bz2
import collections
import contextlib
import gzip
import heapq
import io
import itertools
import operator
import os
import re
import sys
import tempfile

from Bio import SeqIO
import intervaltree
import ncbitax
import zstd

import blastutil


@contextlib.contextmanager
def zstd_open(fname, mode='r'):
    '''Handle both text and byte decompression of the file.'''
    if 'r' in mode:
        with open(fname, 'rb') as fh:
            dctx = zstd.ZstdDecompressor()
            stream_reader = dctx.stream_reader(fh)
            if 'b' not in mode:
                text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
                yield text_stream
                return
            yield stream_reader
    else:
        with open(fname, 'wb') as fh:
            cctx = zstd.ZstdCompressor(level=kwargs.get('level', 10),
                                       threads=kwargs.get('threads', 1))
            stream_writer = cctx.stream_writer(fh)
            if 'b' not in mode:
                text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
                yield text_stream
                return
            yield stream_writer


def compressed_open(fname, mode='r', **kwargs):
    '''Create file-like context manager from an optionally compressed filename.

    Supports regular files, gzip, bz2, and zstd.
    '''

    if fname.endswith('.gz'):
        # Allow using 'level' kwarg as an alias for gzip files.
        if 'level' in kwargs:
            kwargs['compresslevel'] = kwargs.pop('level')
        return gzip.open(fname, mode=mode, **kwargs)
    elif fname.endswith('.bz2'):
        return bz2.open(fname, mode=mode, **kwargs)
    elif fname.endswith('.zst'):
        return zstd_open(fname, mode=mode, **kwargs)
    else:
        return open(fname, mode=mode, **kwargs)



def blast_grouped(m8_file,
                  paired=False,
                  min_bit_score=None,
                  max_expected_value=None):
    min_bit_score = min_bit_score if min_bit_score is not None else 50
    max_expected_value = max_expected_value if max_expected_value is not None else 0.01
    records = blastutil.blast_records(m8_file)
    records = (r for r in records if r.e_val <= max_expected_value)
    records = (r for r in records if r.bit_score >= min_bit_score)
    if paired:
        records = (paired_query_id(rec) for rec in records)
    blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
    return blast_groups


def blast_split_hits(db,
                     m8_file,
                     output=None,
                     paired=False,
                     min_bit_score=None,
                     max_expected_value=None,
                     top_percent=None):
    '''

    Writes tsv output: query_id \t tax_id

    Args:
      db: (TaxonomyDb) Taxonomy db.
      m8_file: (io) Blast m8 file to read.
      output: (io) Output file.
      paired: (bool) Whether to count paired suffixes /1,/2 as one group.
      min_bit_score: (float) Minimum bit score or discard.
      max_expected_value: (float) Maximum e-val or discard.
      top_percent: (float) Only this percent within top hit are used.
    '''
    top_percent = top_percent if top_percent is not None else 10
    all_hits = collections.defaultdict(float)
    for blast_group in blast_grouped(
            m8_file, paired=paired, min_bit_score=min_bit_score, max_expected_value=max_expected_value):
        blast_group = list(blast_group)


        hit_taxids = set()
        hits = [hit for hit in blast_group if hit.taxids]
        if len(hits) == 0:
            continue
        best_score = max(hit.bit_score for hit in hits)
        cutoff_bit_score = (100 - top_percent) / 100 * best_score
        valid_hits = (hit for hit in hits if hit.bit_score >= cutoff_bit_score)
        valid_hits = list(valid_hits)
        # Sort requires realized list
        for hit in valid_hits:
            for txid in hit.taxids:
                hit_taxids.add(txid)

        hit_taxids = list(sorted(list(hit_taxids)))
        n_taxids = len(hit_taxids)
        for i in range(len(hit_taxids)):
            all_hits[hit_taxids[i]] += 1 / n_taxids

        if output:
            query_id = blast_group[0].query_id
            tax_id_s = ':'.join([str(x) for x in hit_taxids])
            print(query_id, tax_id_s, sep='\t', file=output)
    return all_hits


def blast_lca_hits(db,
                   m8_file,
                   output=None,
                   paired=False,
                   min_bit_score=None,
                   max_expected_value=None,
                   top_percent=None):
    '''

    Writes tsv output: query_id \t tax_id

    Args:
      db: (TaxonomyDb) Taxonomy db.
      m8_file: (io) Blast m8 file to read.
      output: (io) Output file.
      paired: (bool) Whether to count paired suffixes /1,/2 as one group.
      min_bit_score: (float) Minimum bit score or discard.
      max_expected_value: (float) Maximum e-val or discard.
      top_percent: (float) Only this percent within top hit are used.
    '''
    top_percent = top_percent if top_percent is not None else 10
    hits = collections.Counter()
    for blast_group in blast_grouped(
            m8_file, paired=paired, min_bit_score=min_bit_score, max_expected_value=max_expected_value):
        blast_group = list(blast_group)
        tax_id = blast_hits_taxid_lca(db, blast_group, top_percent)
        if tax_id is not None:
            hits[tax_id] += 1
        if output:
            query_id = blast_group[0].query_id
            print(query_id, tax_id, sep='\t', file=output)
    return hits


def blast_hits_taxid_lca(db, hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      hits: []BlastRecord groups of hits.
      top_percent: (float) Only consider hits within this percent of top bit score.

    Return:
      (int) Tax id of LCA.
    '''
    hits = [hit for hit in hits if len(hit.taxids)]
    if len(hits) == 0:
        return

    best_score = max(hit.bit_score for hit in hits)
    cutoff_bit_score = (100 - top_percent) / 100 * best_score
    valid_hits = (hit for hit in hits if hit.bit_score >= cutoff_bit_score)
    valid_hits = list(valid_hits)
    # Sort requires realized list
    valid_hits.sort(key=operator.attrgetter('bit_score'), reverse=True)
    if valid_hits:
        tax_ids = list(itertools.chain(*(hit.taxids for hit in valid_hits)))
        return db.coverage_lca(tax_ids, lca_percent=None)


def format_number(num):
    if isinstance(num, int):
        return str(num)
    else:
        return '%.2f' % num

def kraken_dfs(db, lines, taxa_hits, total_hits, taxid, level):
    '''Recursively do DFS for number of hits per taxa.'''
    cum_hits = num_hits = taxa_hits.get(taxid, 0)
    for child_taxid in db.children[taxid]:
        cum_hits += kraken_dfs(db, lines, taxa_hits, total_hits, child_taxid, level + 1)
    percent_covered = cum_hits / total_hits * 100
    rank = ncbitax.kraken_rank_code(db.ranks[taxid])
    name = db.names[taxid]
    if cum_hits > 0:
        lines.append('\t'.join([format_number(percent_covered), format_number(cum_hits), format_number(num_hits), rank,
                                str(taxid), '  ' * level + name]))
    return cum_hits

def kraken_dfs_report(db, taxa_hits, total_reads=None):
    '''Return a kraken compatible DFS report of taxa hits.

    Args:
    db: (TaxonomyDb) Taxonomy db.
    taxa_hits: (collections.Counter) # of hits per tax id.
    total_reads: (int) Total reads (add to unclassified if taxa_hits is less)

    Return:
    []str lines of the report
    '''
    total_hits = sum(taxa_hits.values())
    lines = []
    if total_reads is not None:
        if total_reads >= total_hits:
            unclassified_hits = total_reads - total_hits
            total_hits = total_reads
        elif total_reads < total_hits:
            raise Exception('Total reads given is < number of hits in the report.')
    else:
        unclassified_hits = 0

    unclassified_hits += taxa_hits.get(0, 0)
    unclassified_hits += taxa_hits.get(-1, 0)

    kraken_dfs(db, lines, taxa_hits, total_hits, 1, 0)
    if unclassified_hits > 0:
        percent_covered = unclassified_hits / total_hits * 100
        lines.append(
            '\t'.join([
                format_number(percent_covered), format_number(unclassified_hits), format_number(unclassified_hits), 'U', '0', 'unclassified'
            ])
        )
    return reversed(lines)

def blast_report(args):
    assert args.blast_report or args.blast_taxids
    tax_db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)
    with contextlib.ExitStack() as ctx:
        if args.blast_report and not args.blast_taxids:
            _, blast_taxids_fn = tempfile.mkstemp('.blast.taxids.tsv')
        else:
            blast_taxids_fn = args.blast_taxids

        blast_m8_f = ctx.enter_context(compressed_open(args.blast_m8, 'rt'))
        blast_taxids_f = ctx.enter_context(compressed_open(blast_taxids_fn, 'wt'))
        if args.mode == 'lca':
            hits = blast_lca_hits(tax_db, blast_m8_f, blast_taxids_f, min_bit_score=args.min_bit_score,
                                  max_expected_value=args.max_expected_value, top_percent=args.top_percent)
        elif args.mode == 'split':
            hits = blast_split_hits(tax_db, blast_m8_f, blast_taxids_f, min_bit_score=args.min_bit_score,
                                    max_expected_value=args.max_expected_value, top_percent=args.top_percent)

        if not args.blast_report:
            return

        blast_report_f = ctx.enter_context(compressed_open(args.blast_report, 'wt'))
        if blast_report_f:
            for line in kraken_dfs_report(tax_db, hits, total_reads=args.total_reads):
                print(line, file=blast_report_f)

def add_report_command(subparsers):
    parser = subparsers.add_parser('report')
    parser.add_argument('blast_m8', help='Input BLAST m8 file.')
    parser.add_argument('--blast-report')
    parser.add_argument('--blast-taxids')
    parser.add_argument('--total-reads')
    parser.add_argument('--mode', default='lca', choices=('lca', 'split'), help='Mode to assign hits to taxids')
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--min-bit-score', type=float, default=70)
    parser.add_argument('--max-expected-value', type=float, default=1e-6)
    parser.add_argument('--top-percent', type=float, default=10)
    parser.add_argument('--min-identity', default=0.90, type=float, help='Only consider blast hits with >x identity')
    # parser.add_argument('--min-coverage', default=0.3, type=float, help='Filter out hits with >x coverage ratio')
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=blast_report)

def blast_filter(args):
    # assert args.blast_report or args.blast_lca
    # tax_db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)
    min_bit_score = args.min_bit_score
    max_expected_value = args.max_expected_value
    top_percent = args.top_percent
    with contextlib.ExitStack() as ctx:
        blast_m8_f = ctx.enter_context(compressed_open(args.blast_m8, 'rt'))
        if args.output :
            output_f = ctx.enter_context(compressed_open(args.output, 'wt'))
        else:
            output_f = sys.stdout
        records = blastutil.blast_records(blast_m8_f)
        records = (r for r in records if r.e_val <= max_expected_value)
        records = (r for r in records if r.bit_score >= min_bit_score)
        if args.paired:
            records = (paired_query_id(rec) for rec in records)
        blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
        for blast_group in blast_groups:
            hits = list(blast_group)

            best_score = max(hit.bit_score for hit in hits)
            cutoff_bit_score = (100 - top_percent) / 100 * best_score
            valid_hits = list(hit for hit in hits if hit.bit_score >= cutoff_bit_score and hit.percent_identity >= args.min_identity * 100)
            # Sort requires realized list
            valid_hits.sort(key=operator.attrgetter('bit_score'), reverse=True)
            for hit in valid_hits:
                print(hit, file=output_f)


def add_filter_command(subparsers):
    parser = subparsers.add_parser('filter')
    parser.add_argument('blast_m8', help='Input BLAST m8 file.')
    parser.add_argument('--output')
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--min-bit-score', type=float, default=70)
    parser.add_argument('--max-expected-value', type=float, default=1e-6)
    parser.add_argument('--top-percent', type=float, default=10)
    parser.add_argument('--min-identity', default=0.90, type=float, help='Only consider blast hits with >x identity')
    # parser.add_argument('--min-coverage', default=0.3, type=float, help='Filter out hits with >x coverage ratio')
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=blast_filter)


def blast_truncate(args):
    '''
    Truncates the last query name's blast hits from a blast m8 tab separated
    file. If a previous blast execution was ended before finishing properly, we
    can't trust the hits of the last query and it's best to truncate the file
    off. Use some file seeking tricks to make this much faster than processing
    the whole file.
    '''
    last_qname = None
    fn = args.blast_m8_fn
    if not os.path.exists(fn):
        return
    if os.stat(fn).st_size < 1024:
        with open(fn, 'r+b') as f:
            f.truncate()
        return

    with open(fn, 'r+b') as f:
        f.seek(-1024, 2)
        f.readline()
        # First find the last query name
        for line in f.readlines():
            line = line.decode('ascii')
            parts = line.split('\t')
            qname = parts[0]
            last_qname = qname

        f.seek(0, 2)
        seek_bytes = -2048
        # Seek backwards looking for an earlier query name (!= last query name)
        while True:
            try:
                f.seek(seek_bytes, 1)
            except OSError:  # Hit beginning of file, just truncate to beginning
                f.seek(0)
                f.truncate()
                return

            txt = f.readline()
            if f.readline().decode('ascii').split('\t')[0] != last_qname:
                break

        # We're at an earlier query name, now read lines forward until we find
        # the last query name again. Then un-seek the first line of the last
        # query name and truncate there
        while True:
            line = f.readline()
            if line.decode('ascii').split('\t')[0] == last_qname:
                f.seek(-len(line), 1)
                f.truncate(f.tell())
                print(last_qname)
                return

        # We should never reach here - somehow didn't find the last query name?
        sys.exit(1)


def add_truncate_command(subparsers):
    parser = subparsers.add_parser('truncate')
    parser.add_argument('blast_m8_fn')
    parser.set_defaults(func=blast_truncate)


def blast_fasta_filter(args):
    # max_expected_value = 0.01
    # min_bit_score = 50
    remove_qnames = set()
    with compressed_open(args.blast_m8_fn, 'rt') as f:

        records = blastutil.blast_records(f)
        blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
        for blast_group in blast_groups:
            blast_group = list(blast_group)
            itree = intervaltree.IntervalTree()
            for r in blast_group:
                if r.percent_identity / 100 < args.min_identity:
                    continue

                start = min([r.subject_start, r.subject_end])
                stop = max([r.subject_start, r.subject_end])
                itree[start:stop] = r.percent_identity

            itree.merge_overlaps()
            covered_len = 0
            for iv in itree:
                covered_len += iv.end - iv.begin

            query_len = int(blast_group[0].query_id.split('_')[3])

            covered = covered_len / query_len
            if covered >= args.min_coverage:
                remove_qnames.add(blast_group[0].query_id)

    with open(args.output_fasta, 'wt') as out_f:
        for record in SeqIO.parse(args.fasta, 'fasta'):
            if record.id in remove_qnames:
                continue
            SeqIO.write(record, out_f, 'fasta')


def add_fasta_filter_command(subparsers):
    parser = subparsers.add_parser('fasta-filter')
    parser.add_argument('blast_m8_fn')
    parser.add_argument('-o', '--output-fasta', required=True, help='Output filtered contigs fasta.')
    parser.add_argument('-f', '--fasta', required=True, help='Contigs fasta to filter')
    parser.add_argument('--min-identity', default=0.65, type=float, help='Only consider blast hits with >x identity')
    parser.add_argument('--min-coverage', default=0.3, type=float, help='Filter out hits with >x coverage ratio')
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=blast_fasta_filter)



def spades_key(x):
    return int(x.split('_')[1])


def overlap(start, stop, ivl):
    l = stop - start
    if ivl.begin <= start and ivl.end >= stop:
        return 1
    elif ivl.begin >= start and ivl.end <= stop:
        return (ivl.end - ivl.begin) / l
    elif ivl.begin <= start:
        return (ivl.end - start) / l
    elif ivl.end >= stop:
        return (stop - ivl.begin) / l


def simple_sort(hits):
    # hits.sort(key=lambda r: r.e_val)
    hits.sort(key=lambda r: r[10])
    return hits


def per_species(hits):
    seen_taxids = set()
    for r in hits:
        if all(taxid in seen_taxids for taxid in r[14]):
            continue
        else:
            for taxid in r[14]:
                seen_taxids.add(taxid)
            yield r


def viral_hits1(db, hits):
    # hits.sort(key=lambda r: r.e_val)
    hits.sort(key=lambda r: r[10])
    for r in hits:
        lca_taxid = db.coverage_lca(r.taxids)
        path = db.parent_path(lca_taxid)
        if VIRUSES_TAXID in path and HEV_TAXID not in path:
            yield r


def filter_taxids1(context, hits):
    itree = intervaltree.IntervalTree()
    n_skipped = 0
    output_hits = []
    for r in hits:

        skip_hit = False
        # start = min([r.query_start, r.query_end])
        # stop = max([r.query_start, r.query_end])
        # olap = 0
        # for ivl in itree[start:stop]:
        #     olap = overlap(start, stop, ivl)
        #     if olap >= 0.9 and ivl.data >= r.bit_score:
        #         skip_hit = True
        #         break

        # if skip_hit:
        #     n_skipped += 1
        #     if n_skipped % 1000 == 0:
        #         print('N skipped: {}'.format(n_skipped))
        #     continue
        # itree[start:stop] = r.bit_score


        if not any(taxid in context.whitelist_taxids for taxid in r.taxids):
            continue

        if context.blacklist_taxids:
            for taxid in r[14]:
                if taxid in context.blacklist_taxids:
                    skip_hit = True
                    break
        if skip_hit:
            continue
        output_hits.append(r)
    if output_hits:
        output_hits.sort(key=lambda x: x.bit_score, reverse=True)
        min_bit_score = output_hits[-1].bit_score
        for r in hits:
            if r.bit_score >= min_bit_score - 50:
                yield r


def qgroups_to_hits(qgroups):
    for k, v in qgroups:
        hits = list(itertools.chain(*(list(igrp[1]) for igrp in v)))
        yield hits


def simple_combine(qgroups, out_f):
    for k, v in qgroups:
        hits = list(itertools.chain(*(list(igrp[1]) for igrp in v)))
        hits.sort(key=lambda r: r.e_val)
        for r in hits:
            print(str(r) + '\t' + r.blast_type, file=out_f)

def viral_hits(db, qgroups, out_f):
    for k, v in qgroups:
        hits = list(itertools.chain(*(list(igrp[1]) for igrp in v)))
        hits.sort(key=lambda r: r.e_val)
        for r in hits:
            lca_taxid = db.coverage_lca(r[14])
            path = db.parent_path(lca_taxid)
            if path and VIRUSES_TAXID in path and HEV_TAXID not in path:
                print(str(r) + '\t' + r.blast_type, file=out_f)


class CombineContext:

    def __init__(self, db):
        self.db = db
        self.blacklist_taxids = None
        self.whitelist_taxids = None
        self.path_cache = {}
        self.lca_cache = {}

    def load_lists(self, args):
        # Only used for filtering
        if args.blacklist_taxids:
            blacklist_taxids = set(int(x) for x in args.blacklist_taxids.split(','))
        else:
            blacklist_taxids = set([PRIMATES_TAXID])
        self.blacklist_taxids = set(ncbitax.collect_children(self.db.children, blacklist_taxids))

        if args.whitelist_taxids:
            whitelist_taxids = set(int(x) for x in args.whitelist_taxids.split(','))
            self.whitelist_taxids = set(ncbitax.collect_children(self.db.children, whitelist_taxids))


def blast_combine(args):
    # db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True)
    db = None

    context = CombineContext(db)
    # context.load_lists(args)

    with compressed_open(args.output, 'wt') as out_f:

        group_itrs = []
        for fn, blast_type in zip([args.megablast, args.blastn, args.blastx], ['megablast', 'blastn', 'blastx']):

            if not fn:
                continue
            f = compressed_open(fn, 'rt')
            records = blastutil.blast_tuples(f, blast_type=blast_type)
            blast_groups = itertools.groupby(records, lambda x: x[0])
            group_itrs.append(blast_groups)


        qgroups = itertools.groupby(heapq.merge(*group_itrs, key=lambda x: spades_key(x[0])), lambda x: x[0])

        hits_gen = qgroups_to_hits(qgroups)

        # Hits is generator of hits per each query id
        n_recs = 0
        for i, hits in enumerate(hits_gen):
            n_recs += len(hits)
            if i % 100 == 0:
                print('%s records, %s contigs processed' % (n_recs, i), file=sys.stderr)
            if args.action == 'combine':
                out_hits = hits
            elif args.action == 'viruses':
                out_hits = viral_hits1(db, hits)
            elif args.action == 'filter':
                out_hits = filter_taxids1(context, hits)
            elif args.action == 'per_species':
                out_hits = per_species(hits)
            for r in out_hits:
                r[14] = ';'.join(str(x) for x in r[14]) or 'N/A'
                print('\t'.join(str(x) for x in r), file=out_f)
                # print(str(r).rstrip() + '\t' + r.blast_type, file=out_f)


def add_combine_command(subparsers):
    parser = subparsers.add_parser('combine')
    parser.add_argument('--output')
    parser.add_argument('--megablast', help='Megablast m8')
    parser.add_argument('--blastn', help='Blastn m8')
    parser.add_argument('--blastx', help='Blastx m8')
    parser.add_argument('--blacklist-taxids', help='Comma separated taxids to blacklist - takes precedence over whitelist')
    parser.add_argument('--whitelist-taxids', help='Comma separated taxids to whitelist')
    parser = ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=blast_combine)


def dustmasked(fn):
    with open(fn, 'rt') as f:
        masked = []
        for line in f:
            if line.startswith('>'):
                if masked:
                    yield (qname, masked)
                    masked = []
                qname = line[1:].rstrip()
                continue
            mo = re.match(r'(\d+) - (\d+)', line)
            start = int(mo.group(1))
            end = int(mo.group(2))
            masked.append((start, end))
        yield (qname, masked)


def dust_filter(args):
    dusts = {}
    unmasked_bp = {}
    for qname, masked in dustmasked(args.dust):
        dusts[qname] = masked
        last_end = 0
        # Doesn't count the final unmasked interval
        unmasked_len = 0
        for start, end in masked:
            if start - last_end >= args.min_fragment:
                unmasked_len += start - last_end
            last_end = end
        unmasked_bp[qname] = unmasked_len


    with open(args.output, 'wt') as out_f:
        for record in SeqIO.parse(args.input, 'fasta'):
            seq_len = len(record)
            if seq_len < args.min_len:
                print('Filtered out contig {} - too short'.format(record.id), file=sys.stderr)
                continue
            dust = dusts.get(record.id)
            if dust:
                last_ivl = dust[-1]
                total_len = unmasked_bp[record.id]
                last_len = seq_len - last_ivl[1]
                if last_len >= args.min_fragment:
                    total_len += last_len
                lc_prop = total_len / seq_len
                if lc_prop >= 0.2:
                    print('Kept contig {} - {:2f}% unique content'.format(record.id, lc_prop * 100), file=sys.stderr)
                else:
                    print('Filtered contig {} - {:2f}% unique content'.format(record.id, lc_prop * 100), file=sys.stderr)
                    continue

            SeqIO.write(record, out_f, 'fasta')


def add_dust_filter_command(subparsers):
    parser = subparsers.add_parser('dust-filter')
    parser.add_argument('dust')
    parser.add_argument('-i', '--input', help='Input fasta file')
    parser.add_argument('-o', '--output', help='Filtered output fasta file')
    parser.add_argument('--min-len', help='Minimum length', type=int, default=101)
    parser.add_argument('--min-fragment', help='Minimum fragment length to count between masked intervals', type=int, default=101)
    parser.set_defaults(func=dust_filter)


def main():
    parser = argparse.ArgumentParser(description='Blast toolkit')
    subparsers = parser.add_subparsers()
    add_truncate_command(subparsers)
    add_fasta_filter_command(subparsers)
    add_filter_command(subparsers)
    add_report_command(subparsers)
    add_combine_command(subparsers)
    add_dust_filter_command(subparsers)
    args = parser.parse_args()
    func = getattr(args, 'func', None)
    if func:
        args.func(args)
    else:
        print('No command selected')


if __name__ == '__main__':
    main()
