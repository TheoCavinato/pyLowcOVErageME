from pysam import VariantFile
import argparse

"""
When using high coverage, the PL of a sample should always correspond to the
at a certain SNP should correspond to its corresponding genotype
"""

#=====================
# User's parameters
#=====================

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
args = parser.parse_args()

#=====================
# Read vcf
#=====================

vcf_path = args.vcf
vcf_reader = VariantFile(vcf_path, 'r')
samples = list(vcf_reader.header.samples)
samples_error_PL_GT = [0 for _ in samples]
samples_available_snps = [0 for _ in samples]
for rec in vcf_reader:
    for n, sample in enumerate(samples):
        if None not in rec.samples[sample]['GT']:
            samples_available_snps[n]+=1
            GT = sum(rec.samples[sample]['GT'])
            PL = rec.samples[sample]['PL']
            PL_to_GENO = PL.index(min(PL))
            if PL_to_GENO!=GT:
                print(PL, GT)
                samples_error_PL_GT[n]+=1

vcf_reader.close()

#=====================
# Output number of errors
# for each sample
#=====================
for n in range(len(samples)):
    print(samples[n], samples_error_PL_GT[n]/samples_available_snps[n])