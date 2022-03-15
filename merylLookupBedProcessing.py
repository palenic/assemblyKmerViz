import argparse
import merylLookupSuite as mer

parser = argparse.ArgumentParser(description="""Takes a bed file produced by merylLookupWigToBed
    and turns it into a .csv file containing for
    each contig its name, length, and counts of how many
    positions in the contig have that multiplicity,
    optionally scaled based on contig length.
    Scores >4 are aggregated as 5.
    Produces a .csv file with columns: contigName, contigLength, 
    counts_0, counts_1, counts_2, counts_3, counts_4, counts_5_or_more.
    The file is sorted by contigName.""",
    epilog="""Author: Chiara Paleni""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_bed", help="path to input bed file")
parser.add_argument("output_csv", help="path to output csv file")
parser.add_argument("-r", "--relative-counts",help="return the percentage of bases with that multiplicity instead of the counts",
                    action="store_true")
parser.add_argument("-q", "--quiet",help="suppress verbose output",action="store_true")
args = parser.parse_args()
mer.merylLookupBedProcessing(args.input_bed,args.output_csv,args.relative_counts,not args.quiet)

