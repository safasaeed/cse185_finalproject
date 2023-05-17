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
  

genotypes = '


