from pysam import VariantFile
import argparse, numpy

"""
The number of available should respect the coverage we used
"""

#=====================
# User's parameters
#=====================

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--coverage", required=True, type=float)
args = parser.parse_args()

#=====================
# Read vcf
#=====================

vcf_path = args.vcf
vcf_reader = VariantFile(vcf_path, 'r')
samples = list(vcf_reader.header.samples)
samples_sites = [0 for _ in samples]
SNPs=0
for rec in vcf_reader:
    for n, sample in enumerate(samples):
        #=====================
        # Count number of sites
        # available
        #=====================
        if None not in rec.samples[sample]['GT']:
            samples_sites[n] += 1
    SNPs+=1
vcf_reader.close()

expectation=sum([sum([1 for _ in range(SNPs) if numpy.random.poisson(args.coverage, 1)[0]!=0]) for _ in range(len(samples))])/len(samples)
average=sum(samples_sites)/len(samples_sites)

print("##############################################")
print("Sample\tNumberOfSitesAvailable\tNumberOfSitesExpected")
for n in range(len(samples)):
    print(samples[n], '\t', samples_sites[n], '\t', expectation)

print("\n##############################################")
print("Average:", average, "Expectation:", expectation, "Total number of SNPs:", SNPs, '\n')