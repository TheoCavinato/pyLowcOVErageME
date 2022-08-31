from pysam import VariantFile
import argparse, numpy, math

"""
Input a VCF
The algorithm will simulate low coverage data for each individual
Output a new VCF
"""


#=====================
# User's parameters
#=====================

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--coverage", required=True, type=float)
parser.add_argument("--error", required=False, type=float, default=0.0001)
parser.add_argument("--samples", required=False, nargs='+')
args = parser.parse_args()

#=====================
# Useful functions
#=====================

def GL_to_PL(GLs):
    log_values = [-10*GL for GL in GLs]
    min_value = min(log_values)
    PLs = [str(int(round(log_value - min_value))) for log_value in log_values]
    return ','.join(PLs)

def emission_probas(H1, H2, read_allele, error):
    if H1 == H2 == read_allele:
        return math.log10(1 - error)
    elif H1 != read_allele and H2 != read_allele:
        return math.log10(error/3)
    else:
        return math.log10(1/2 - error/3)

def simulate_DP_PL(GT, coverage, error):
    
    # Simulate a number of reads covering the position
    DP = numpy.random.poisson(coverage, 1)[0]
    if DP==0:
        PL = "0,0,0"
        GT_str = "./."

    else:
        GT_str = GT
        GT_int = [int(i) for i in GT.split('|')] if '|' in GT else [int(i) for i in GT.split('/')]

        H0_nbr = numpy.random.binomial(DP, 0.5) #X alleles for hap0
        H1_nbr = DP - H0_nbr #DP - X alleles for hap1

        # Calculate the likelihoods of the genotypes 00, 01 and 11
        p_D_00 = emission_probas( 0, 0, GT_int[0], error)*H0_nbr + emission_probas( 0, 0, GT_int[1], error)*H1_nbr
        p_D_01 = emission_probas( 0, 1, GT_int[0], error)*H0_nbr + emission_probas( 0, 1, GT_int[1], error)*H1_nbr
        p_D_11 = emission_probas( 1, 1, GT_int[0], error)*H0_nbr + emission_probas( 1, 1, GT_int[1], error)*H1_nbr

        PL = GL_to_PL((p_D_00, p_D_01, p_D_11))

    return GT_str + ':' + str(DP)+':'+PL

def output_header(p_vcf_reader):
    # Add PL and DP fields to the header
    split_header = p_vcf_reader.header.__str__().split('\n')
    format_checker = False
    for line in split_header[:-2]:
        if line[:8] == "##FORMAT":
            format_checker=True
        if line[:8] != "##FORMAT" and format_checker:
            print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
            print('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt">')
            format_checker=False
        print(line)
    if format_checker:
            print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
            print('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt">')
    last_line = split_header[-2].split('\t')[:9]
    for sample in list(vcf_reader.header.samples):
        last_line.append(sample)
    print('\t'.join(last_line))

#=====================
# 1. Import the vcf
#=====================

vcf_path = args.vcf
vcf_reader = VariantFile(vcf_path, 'r')
if args.samples:
    vcf_reader.subset_samples(args.samples)

output_header(vcf_reader)

for rec in vcf_reader:

    split_rec = rec.__str__().split()
    split_rec[8]='GT:DP:PL'
    for n, sample_GT in enumerate(split_rec[9:]):

        #=====================
        # 2. simulate low cov
        # for each sample
        #=====================
        split_rec[9+n] = simulate_DP_PL(sample_GT, args.coverage, args.error)

    # Output new line
    print('\t'.join(split_rec))

vcf_reader.close()