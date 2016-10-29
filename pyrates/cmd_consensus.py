"""Commandline interface for consensus calling.
"""

import argparse
import datetime
import resource
import time
import logging

import pyrates.consensus as cons
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
        '--adapter', '-a',
        default='GACT',
        help='Constant part of barcode adapter. This is expected to be located' +
        ' between the UID and the actual read sequence.'
    )
    parser.add_argument(
        '--merge-size', '-m',
        default=3,
        help='Attempt to merge clusters no larger than this into larger clusters.'
    )
    parser.add_argument(
        '--merge-target',
        default=10,
        help='Maximum size for clusters to be considered as a target for the merging' +
        ' of smaller clusters.'
    )
    parser.add_argument(
        '--id-tolerance',
        default=5,
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
    logger.info('This is pyrates ' + __version__)
    logger.debug('At revision ' + get_versions()['full-revisionid'])
    logger.info('Processing input file ' + args.fastq)
    if args.merge_size and args.merge_target:
        logger.info("Will merge clusters smaller than %d into clusters no larger than than %d",
                    args.merge_size, args.merge_target)
    else:
        logger.warning("Merging of small clusters is disabled.")
    logger.info('Consensus sequences will go to ' + args.output)

    output_fun = utils.smart_open(args.output)

    ## start consensus computation
    started_at = time.time()
    seq = cons.from_fastq(args.fastq, args.id_length, args.adapter)

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
        cons_time = time.time()
        logger.info('Time taken for consensus: %s',
                    str(datetime.timedelta(seconds=cons_time - started_at)))

    with output_fun(args.output, 'w') as output:
        if args.merge_size and args.merge_target:
            logger.info('Merging small clusters')
            candidates = []
            targets = []
            merge_count = 0
            for (seq_count, uid) in enumerate(seq):
                if logger.isEnabledFor(logging.DEBUG) and seq_count % 10000 == 0:
                    logger.debug("clusters: %d, merged: %d, small: %d, targets: %d",
                                 seq_count, merge_count, len(candidates), len(targets))
                if seq[uid].size <= args.merge_size:
                    merged = False
                    for consensus in targets:
                        merged = consensus.merge(seq[uid], args.id_tolerance)
                        if merged:
                            merge_count += 1
                            break
                    if not merged:
                        ## attempt merging with other candidates
                        remove = None
                        for (i, cand) in enumerate(candidates):
                            merged = seq[uid].merge(cand, args.id_tolerance)
                            if merged:
                                merge_count += 1
                                remove = i
                                break
                        if merged:
                            del candidates[remove]
                            if seq[uid].size > args.merge_size:
                                targets.append(seq[uid])
                            else:
                                candidates.append(seq[uid])
                        else:
                            candidates.append(seq[uid])
                elif seq[uid].size <= args.merge_target:
                    targets.append(seq[uid])
                    processed = []
                    for (i, cand) in enumerate(candidates):
                        if seq[uid].merge(cand, args.id_tolerance):
                            processed.insert(0, i)
                    for i in processed:
                        del candidates[i]
                else:
                    output.write(str(seq[uid]) + "\n")
            logger.info('Clusters merged: %d', merge_count)
            logger.info('Small clusters remaining: %d', len(candidates))
            for consensus in targets:
                output.write(str(consensus) + "\n")
            for consensus in candidates:
                output.write(str(consensus) + "\n")
            logger.info('Time taken for merging: %s',
                        str(datetime.timedelta(seconds=time.time() - cons_time)))
        else:
            for uid in seq:
                output.write(str(seq[uid]) + "\n")
    logger.info('Total time taken: %s', str(datetime.timedelta(seconds=time.time() - started_at)))
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    logger.info('Memory used: %.2f MB', mem)
