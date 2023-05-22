# cse185_finalproject

This is a demonstration project for CSE185. It implements the linear regression subset of the "plink" package used to perform GWAS (Genome-Wide Association Study). See the official website for plink here: (https://www.cog-genomics.org/plink/)

# Installation Instructions 

Installation requires the `qqman` library to be installed in order to visualize plots. You can install these with `pip`:

```
pip install --user qqman 
```
Once required libraries are installed, you can install `mypileup` with the following command:

```
python setup.py install
```



# Basic usage

The basic usage of `myplink` is:

```
//insert command line 
```

To run `myplink` on a small test example (using files in this repo):
```
//insert test example line
```

This should produce the output below:
```
//insert output 
```

To compare to output of `plink`, run:
```
//compare to lab3 and our own data 
```

# myplink options

The required input to `myplink` is a VCF file (containing the SNPs for a subset of samples) and a .phen file (contains normalized melanin(?) values for each sample) 

* `--pheno FILE`, `--vcf FILE`: //insert info on each 

* `-linear`, //insert blur 

* `-o FILE`, `--output FILE`: Write output to file. By default, output is written to stdout.


# File format

The output file format is a data table as a text file. 

# Contributors

This repository was generated by Safa Saeed, Saara Kriplani, and Honieh Hemati with inspiration from the CSE 185 demo project:  https://github.com/gymreklab/cse185-demo-project/blob/main/mypileup/mypileup.py 

Please submit a pull request with any corrections or suggestions.

# Testing

To run tests:
```
# Run command line tests
//insert 

# Run unit tests
//insert 
```



