"""Commandline interface for consensus calling.
"""

import argparse
import datetime
import resource
import time
import logging

import pyrates.clustering as clust
import pyrates.utils as utils
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
        '--prefix-length', '-p',
        metavar='PREFIX',
        default=5, type=int,
        help="Length of UID prefix to use in read clustering. Larger values may speed up" +
        " the clustering but will require more memory."
    )
    parser.add_argument(
        '--adapter', '-a',
        default='GACT',
        help='Constant part of barcode adapter. This is expected to be located' +
        ' between the UID and the actual read sequence.'
    )
    parser.add_argument(
        '--id-tolerance', '-t',
        default=5, type=int,
        help='Maximum number of differences between IDs allowed to consider merging of clusters.'
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

    ## Add parameter values and version info to log
    logger.info('This is pyrates %s', __version__)
    logger.debug('At revision %s', get_versions()['full-revisionid'])
    logger.info('Processing input file %r', args.fastq)
    logger.info('Consensus sequences will go to %r', args.output)
    logger.info('UID length: %d', args.id_length)
    logger.info('Maximum number of mismatches in UID allowed within cluster: %d', args.id_tolerance)
    logger.info('Length of UID prefix used in clustering: %d', args.prefix_length)
    if args.prefix_length <= args.id_tolerance:
        logger.info('To reduce running time choose a prefix longer than the allowed number' +
                    ' of UID mismatches')
    logger.info('Adapter sequence: %r', args.adapter)


    ## start consensus computation
    started_at = time.time()
    seq = clust.Clustering.from_fastq(input_file=args.fastq, id_length=args.id_length,
                                      adapter=args.adapter, threshold=args.id_tolerance,
                                      prefix=args.prefix_length)
    seq.write(args.output_file)
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
        logger.info("Number of sequences with ambiguous label %d (%.2f%%)",
                    seq.fail_count, seq.fail_count/len(seq)*100)
    logger.info('Total time taken: %s', str(datetime.timedelta(seconds=time.time() - started_at)))
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    logger.info('Memory used: %.2f MB', mem)
