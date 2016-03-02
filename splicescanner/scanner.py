
from maxentpy import maxent
import argparse
import pysam
import os
from collections import namedtuple

Region = namedtuple('Region', ['chr', 'start', 'end'])

def read_regions(region):
    if os.path.exists(region):
        for line in open(region):
            if len(line)==0 or line[0]=='#':
                continue
            toks = line.split()
            yield Region(chr=toks[0], start=int(toks[1]), end=int(toks[2]))
    else:
        toks = region.replace(",", "").split(":")
        if len(toks)<2:
            raise ValueError('Error parsing input region (not a path, and no :')
        startend = toks[1].split('-')
        if len(startend)==2:
            start = int(startend[0])
            end = int(startend[1])
        else:
            start = int(startend[0])
            end = start + 1
        yield Region(chr=toks[0], start=start, end=end)


def emit_scores(ref, region):
    scorer = maxent.SpliceScorer()
    seq_offset = max(maxent.DONOR_JUNCTION_OFFSET, maxent.ACCEPTOR_JUNCTION_OFFSET)
    seq = ref.fetch(region.chr, region.start-seq_offset, region.end + maxent.ACCEPTOR_BASES)
    for offset in range(region.end - region.start):
        donor_score = scorer.score5(seq[offset-maxent.DONOR_JUNCTION_OFFSET+seq_offset:offset+maxent.DONOR_BASES-maxent.DONOR_JUNCTION_OFFSET+seq_offset])
        acceptor_score = scorer.score3(seq[offset-maxent.ACCEPTOR_JUNCTION_OFFSET+seq_offset:offset+maxent.ACCEPTOR_BASES-maxent.ACCEPTOR_JUNCTION_OFFSET+seq_offset])
        print "\t".join([region.chr, str(region.start+offset), str(donor_score), str(acceptor_score)])

def main(args):
    ref = pysam.FastaFile(args.fasta)
    for region in read_regions(args.regions):
        emit_scores(ref, region)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Reference fasta")
    parser.add_argument("-r", "--regions", help="Genomic region chr:start-end or bed file")

    args = parser.parse_args()
    main(args)

