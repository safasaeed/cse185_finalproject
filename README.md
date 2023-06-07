# cse185_finalproject

This is a demonstration project for CSE185. It implements the linear regression subset of the "plink" package used to perform GWAS (Genome-Wide Association Study). See the official website for plink here: (https://www.cog-genomics.org/plink/)

# Installation Instructions 

Installation requires the `qqman` library to be installed in order to visualize plots. You can install these with `pip`:

```
pip install --user qqman 
```
Installation also requires 'PyVCF', which can be installed following the below commands:

```
pip install setuptools==58
pip install pyvcf
```

Please download the setup.py file and the mygwas directory from this Github page. Once required libraries and files are installed, you can install `mypileup` with the following command:

```
python setup.py install
```

You may have to install other dependencies such as numpy, pandas, and statsmodels once running the mygwas command. This can be done using `pip install`.

# Basic usage

The basic usage of `mygwas` is:

```
mygwas [-h] [-o FILE] vcf phen 
```

To run `mygwas` on a small test example (using files in this repo):
```
mygwas simple.vcf simple.phen
```

This should produce a .assoc.linear file as output, which contains chromosome information, SNP names, and p-values.

To compare to output of `plink`, run:
```
plink --pheno FILE --allow-no-sex --linear --out FILE --vcf FILE
```

# mygwas options

The required input to `mygwas` is a VCF file (containing the SNPs for a subset of samples) and a .phen file (contains normalized phenotype values for each sample) 

* `pheno FILE`, `vcf FILE`: Input phenotype and VCF files 

* `-o FILE`, `--output FILE`: Write output to file. By default, output is written to stdout.


# File format

There are two output files, an .assoc.linear file which contains the results of the linear regression and a plot named QQMan.png. 

# Contributors

This repository was generated by Safa Saeed, Saara Kriplani, and Honieh Hemati with inspiration from the CSE 185 demo project:  https://github.com/gymreklab/cse185-demo-project/blob/main/mypileup/mypileup.py 

Please submit a pull request with any corrections or suggestions.

# Testing (to be updated)

To run tests:
```
# Run command line tests
//insert 

# Run unit tests
//insert 
```



