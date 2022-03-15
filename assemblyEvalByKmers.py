# call meryl directly, do not save intermediate .wig output
import argparse
import merylLookupSuite as mer

parser = argparse.ArgumentParser(description="""Evaluate a genome assembly by producing a 
    .bed file containing the multiplicity of k-mers at each position.
    Calls meryl-lookup wig-count 
Summarizes variableStep .wig files produced by 
    meryl-lookup -wig-count showing the multiplicity of the kmer
    starting at each position in the genome assembly. Produces a tab-separated 
    .bed file. The file is 0-based, half-open [start-1,end) and has fields:
    chrom	start	end .	multiplicity    .   start   start   color
    Default:
        Produces a tab-separated .bed file that aggregates adjacent positions with 
        the same multiplicity. Scores >4 are aggregated as 5.
    Binned mode (-b):
        Computes multiplicity inside bins of the specified width (default -b 100)
        with majority vote of all positions inside the bin.
    """,
    epilog="""Author: Chiara Paleni""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_wig", help="path to input wig file",nargs='?', type=argparse.FileType('r'))
parser.add_argument("output_bed", help="path to output bed file")
parser.add_argument("-b","--binned", help="turn on binned mode",action="store_true")
parser.add_argument("-w", help="bin width for binned mode",nargs="?",type=int,default=100)
parser.add_argument("-q", "--quiet",help="suppress verbose output",action="store_true")
args = parser.parse_args()
if(args.binned):
    mer.merylLookupWigToBedBinned(args.input_wig,args.output_bed, args.w,not args.quiet)
else:
    mer.merylLookupWigToBed(args.input_wig,args.output_bed,not args.quiet)

