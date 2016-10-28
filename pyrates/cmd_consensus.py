"""Commandline interface for consensus calling.
"""

import argparse
import datetime
import gzip
import resource
import time
import logging

import pyrates.consensus as cons
import pyrates.utils as utils
import pyrates.sequence as pseq
from . import __version__
from ._version import get_versions

def main():
    """Entrypoint for command-line interface
    """
    ## parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Correct errors in UID labelled reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('fastq',
                        help='Fastq file with input reads')
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        required=True,
        help='Output file name'
        )
    parser.add_argument(
        '--id-length', '-b',
        metavar='LENGTH',
        default=8,
        help='Length of unique identifier at start of read.'
    )
    parser.add_argument(
        '--adapter', '-a',
        default='GACT',
        help='Constant part of barcode adapter. This is expected to be located' +
        ' between the UID and the actual read sequence.'
    )
    parser.add_argument(
        '--log',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        type=str.upper,
        help='Set verbosity of logging output.'
        )
    parser.add_argument(
        '--version', '-V', action='version',
        version='%(prog)s ' + __version__
    )
    args = parser.parse_args()

    ## configure logging
    logger = utils.get_logger('pyrates', args.log, [utils.console_handler()])
    logger.info('This is pyrates ' + __version__)
    logger.debug('At revision ' + get_versions()['full-revisionid'])
    logger.info('Processing input file ' + args.fastq)
    logger.info('Consensus sequences will go to ' + args.output)

    input_fun = open
    if args.fastq.endswith('.gz'):
        input_fun = gzip.open
    output_fun = open
    if args.output.endswith('.gz'):
        output_fun = gzip.open

    ## start consensus computation
    started_at = time.time()
    line_count = 0
    total_skipped = 0
    seq = {}
    id_length = args.id_length
    adapt_length = args.id_length + len(args.adapter)
    with input_fun(args.fastq) as fastq:
        for line in fastq:
            # print out some stats as we go
            if logger.isEnabledFor(logging.DEBUG) and (line_count % 100000) == 0:
                logger.debug("reads: %d clusters: %d skipped: %d",
                             line_count/4, len(seq), total_skipped)
            elif (line_count % 4) == 1:
                line = line.rstrip("\n")
                nameid = line[0:id_length] + line[-id_length:]
                sequence = line[adapt_length:-adapt_length]
            elif (line_count % 4) == 3:
                line = line.rstrip("\n")
                qnameid = line[0:id_length] + line[-id_length:]
                qsequence = line[adapt_length:-adapt_length]

                uid = pseq.SequenceWithQuality(nameid, qnameid)
                read_seq = pseq.SequenceWithQuality(sequence, qsequence)
                if nameid not in seq:
                    seq[nameid] = cons.Consensus(uid, read_seq)
                else:
                    success = seq[nameid].update(uid, read_seq)
                    if not success:
                        total_skipped += 1
            line_count += 1

    if logger.isEnabledFor(logging.INFO):
        total_different = 0
        total_shorter = 0
        total_longer = 0

        for uid in seq:
            total_different += seq[uid].different
            total_shorter += seq[uid].shorter
            total_longer += seq[uid].longer
        logger.info("Number of consensus sequence with unique labels: %d", len(seq))
        logger.info("Number sequences grossly differet from consensus with same label: %d",
                    total_different)
        logger.info("Number of sequences that were shorter than consensus sequence: %d",
                    total_shorter)
        logger.info("Number of sequences that were longer then consensus sequence %d",
                    total_longer)

    #
    # Print everything out
    # '@'int [int,A,int,C,int,T,int,G,int,N,int ...]
    #
    with output_fun(args.output, 'w') as output:
        for consensus in sorted(seq.values(), key=lambda x: x.size):
            output.write(str(consensus) + "\n")

    logger.info('Time taken: %s', str(datetime.timedelta(seconds=time.time() - started_at)))
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    logger.info('Memory used: %.2f MB', mem)
