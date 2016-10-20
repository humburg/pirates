"""Commandline interface for consensus calling.
"""

import argparse
import datetime
import gzip
import logging
import resource
import time

import pyrates.consensus as cons

LOG = logging.getLogger(__name__)

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
        '--log',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        type=str.upper,
        help='Set verbosity of log file output.'
        )
    args = parser.parse_args()

    ## configure logging
    LOG.setLevel(getattr(logging, args.log, None))
    console_handler = logging.StreamHandler()
    cons_formatter = logging.Formatter('[%(levelname)s] %(name)s - %(message)s')
    console_handler.setFormatter(cons_formatter)
    LOG.addHandler(console_handler)

    input_fun = open
    if args.fastq.endswith('.gz'):
        input_fun = gzip.open
    output_fun = open
    if args.output.endswith('.gz'):
        output_fun = gzip.open

    ## start consensus computation
    started_at = time.time()
    line_count = 0
    seq = {}
    seqcount = {}
    different = {}
    shorter = {}
    longer = {}
    diff = {}
    with input_fun(args.fastq) as fastq:
        for line in fastq:
            # print out some stats as we go
            if (line_count % 100000) == 0:
                LOG.debug(line_count, "seen", len(seq),
                          "different", len(different), "shorter", len(shorter),
                          "longer", len(longer))

            if (line_count % 4) == 1:
                line = line.rstrip("\n")
                nameid = line[0:8] + line[-8:]
                sequence = line[12:-12]
                line_count += 1
                continue

            if (line_count % 4) == 3:
                line = line.rstrip("\n")
                qnameid = line[0:8] + line[-8:]
                qsequence = line[12:-12]
                line_count += 1

                # if not already seen then record and move on, create zero diffs dict
                if nameid not in seq:
                    seq[nameid] = [qnameid, sequence, qsequence]
                    seqcount[nameid] = 1
                    diff[nameid] = {}
                    continue

                # we have this id already recorded so must process for possible consensus
                #
                # if grossly different then just count this and move on
                if cons.grosslydifferent(sequence, seq[nameid][1]):
                    if nameid not in different:
                        different[nameid] = 0

                    different[nameid] += 1
                    continue

                # if sequence length is shorter, count this occurance, abandon this
                # sequence and move on
                if len(sequence) < len(seq[nameid][1]):
                    if nameid not in shorter:
                        shorter[nameid] = 0

                    shorter[nameid] += 1
                    continue

                # if new sequence is longer, count this occurance
                # replace consensus sequence if built from only one other sequence
                if len(sequence) > len(seq[nameid][1]):
                    if nameid not in longer:
                        longer[nameid] = 0
                    longer[nameid] += 1

                    # if consensus built from just one sequence then replace
                    if seqcount[nameid] == 1:
                        seq[nameid] = [qnameid, sequence, qsequence]
                        continue
                    # else abandon this sequence
                    continue

                # do consensus
                #
                result = cons.consensus(seq[nameid][0], qnameid, seq[nameid][1],
                                        sequence, seq[nameid][2], qsequence,
                                        seqcount[nameid], diff[nameid])

                # {"qid": string, "seq": string, "qseq": string}
                seq[nameid][0] = result['qid']
                seq[nameid][1] = result['seq']
                seq[nameid][2] = result['qseq']

                # update count for this sequence label
                seqcount[nameid] += 1
                continue

            line_count += 1

    LOG.info("Number of consensus sequence with unique labels " + str(len(seq)))
    LOG.info("Number sequences grossly difference from consensus with same label " +
             str(len(different)))
    LOG.info("Number of sequences that were shorter than consensus sequence " +
             str(len(shorter)))
    LOG.info("Number of sequences that were longer then consensus sequence " + str(len(longer)))

    #
    # Print everything out
    # '@'int [int,A,int,C,int,T,int,G,int,N,int ...]
    #
    with output_fun(args.output, 'w') as output:
        for label in sorted(seqcount, key=lambda x: seqcount[x]):
            name = "@" + str(seqcount[label])
            # loop over any diffs and add to name string
            for pos in diff[label]:
                name += " " + str(pos)
                for nuc in ['A', 'C', 'T', 'G', 'N']:
                    name += nuc + str(diff[label][pos][nuc])

            # print out sequence data
            output.write(name + "\n")
            output.write(label + seq[label][1] + "\n+\n")
            output.write(seq[label][0] + seq[label][2] + "\n")

    LOG.info('Time taken: ' + str(datetime.timedelta(seconds=time.time() - started_at)))
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    LOG.info('Memory used: ' + str(mem) + 'MB')
