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
    parser.add_argument("vcf", help = "Input unzipped vcf file with genotype data", type = str)
    # Phenotype input file
    parser.add_argument("phen", help = "Input phen file with phenotype data", type = str)

    # Output
    parser.add_argument("-o", "--out", help = "Specify output file name", metavar = "FILE", type = str, required = False)

    # Parse args
    args = parser.parse_args()

    # Open VCF file
    genotypes = args.vcf
    vcf_reader = vcf.Reader(open(genotypes, 'r'))
    
    # Open phen file
    phenotypes = pd.read_csv(args.phen, sep = "\t", header = None)

    # Access header info
    header = vcf_reader.metadata
    sample_names = vcf_reader.samples

    # Iterate over SNPs in VCF file
    for record in vcf_reader:
        sample_data = pd.DataFrame(columns = ['Sample', 'Genotype', 'Phenotype'])
        for call in record.samples:
            sample_name = call.sample
            gt = call.data.GT
            allele1, allele2 = gt.split("|")
            gt_val = int(allele1) + int(allele2)
            pt = phenotypes.loc[phenotypes[0] == sample_name, 2].item()
            sample_data = sample_data.append({'Sample' : sample_name, 'Genotype' : gt_val, 'Phenotype' : pt}, ignore_index=True)
   
    # Convert Genotype and Phenotype columns to numeric type
    sample_data['Genotype'] = pd.to_numeric(sample_data['Genotype'])
    sample_data['Phenotype'] = pd.to_numeric(sample_data['Phenotype'])
    
    # Perform linear regression between phenotype (y) and genotype (x)
    x = sample_data['Genotype']
    y = sample_data['Phenotype']
    
    x = sm.add_constant(x)
    model = sm.OLS(y,x)
    results = model.fit()
    
    # Print GWAS results
    print(results.summary())
