# pyLowcOVErageME
Convert a VCF to a low coverage VCF with the same samples.
## Quick start
For instance, if you have a vcf and you want it to be 0.1x coverage you can run:\
`python3 sim_low_cov.py --vcf your_vcf.vcf.gz --coverage 0.1` \
If you want to directly compress the outptut just do: \
`python3 sim_low_cov.py --vcf your_vcf.vcf.gz --coverage 0.1 | bgzip -c > ouptut.vcf.gz`

## Options
`--samples` subsample the vcf to a specified list of samples \
Example: `python3 sim_low_cov.py --vcf you_vcf.vcf.gz --coverage 0.1 --samples NA20792 NA20795 NA20796` \
`--error` change the error rate used in the PL calculation. Default is 0.0001