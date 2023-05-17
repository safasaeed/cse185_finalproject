import numpy as np
import pandas as pd
import statsmodels.api as sm
import vcf
import argparse

def main():
    parser = argparse.ArgumentParser(
        prog = "mygwas",
        description = "Command-line script to perform GWAS linear regression"
    )

    # Input
    # Genotype input file
    parser.add_argument("vcf", help = "Input vcf file with genotype data", type = str)
    # Phenotype input file
    parser.add_argument("phen", help = "Input phen file with phenotype data", type = str)

    # Output
    parser.add_argument("-o", "--out", help = "Specify output file name", metavar = "FILE" type = str, required = False)

    # Parse args
    args = parser.parse_args()

    # Open VCF file
    genotypes = args.vcf
    vcf_reader = vcf.Reader(open(genotypes, 'r'))
    
    # Open phen file
    phenotypes = args.phen
    

    # Access header info
    header = vcf_reader.metadata
    sample_names = vcf_reader.samples

    # Iterate over SNPs in VCF file
    for record in vcf_reader:
        sample_data = []
        for call in record.samples:
            sample_name = call.sample
            gt = call.data.GT
