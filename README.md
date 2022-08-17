# pyLowcOVErageME
Convert a VCF to a low coverage VCF with the same samples.
You have a vcf with and you want it to be 0.1x coverage? Do
`python3 sim_low_cov.py --vcf your_vcf.vcf.gz --coverage 0.1`
If you want to directly compress the outptut just do:
`python3 sim_low_cov.py --vcf your_vcf.vcf.gz --coverage 0.1 | bgzip -c > ouptut.vcf.gz`