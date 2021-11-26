KPBS
====

kpbs is a python3 script for sliding-window analysis of population branch statistics (PBS). The sliding-window in kpbs is based on nucleotide position in chromosomes/scaffolds not number of SNP positions.

## Requirement

kpbs need following external libraries: [NumPy](https://pypi.org/project/numpy/), [scikit-allel](https://pypi.org/project/scikit-allel/), [pandas](https://pypi.org/project/pandas/)

### Recommended operating environment

- Python ≥ 3.6
- NumPy ≥ 1.19
- scikit-allel ≥ 1.3
- pandas ≥ 1.3

## Usage

### Simple sliding-window PBS analysis

Specify VCF file, three files for lists of IDs of each population and parameters for sliding-window (window size and step size) 

```
python3 kpbs.py --vcf [VCF file] --pop1 [file for list of pop1 IDs] --pop2 [file for list of pop2 IDs] --pop3 [file for list of pop3 IDs] -w [window size (integer)] -s [step size (integer)]
```

### Sliding-window PBS analysis with *P*-value calculations using intergenic regions (needs GFF3 file)

Specify number for bootstrap iteration (more than 1,000 is recommended) and GFF3 file in addition to the simple analysis above

```
python3 kpbs.py --vcf [VCF file] --pop1 [file for list of pop1 IDs] --pop2 [file for list of pop2 IDs] --pop3 [file for list of pop3 IDs] -w [window size (integer)] -s [step size (integer)] --num_bs [number of bootstrap iteration (integer)] --gff [GFF3 file] 
```

### Arguments

```
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     vcf file
  -p1 POP1, --pop1 POP1
                        sample ID list for population1
  -p2 POP2, --pop2 POP2
                        sample ID list for population2
  -p3 POP3, --pop3 POP3
                        sample ID list for population3
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        window size in nucleotide length (bp) to be used in the window analysis 
  -s STEP_SIZE, --step_size STEP_SIZE
                        step size in nucleotide length (bp) to be used in the window analysis
  --gff GFF             GFF3 file
  --num_bs NUM_BS       number of bootstrap iteration (integer)
```



## Licence

The source code is licensed [MIT](https://github.com/takuronkym/kpbs/blob/58516b92c593ea52892b35803b61227fa0040a15/LICENSE).

## Authors

This script was developed by Takuro Nakayama and Atsushi Ikemoto from [laboratory for evolutionary biology](https://klabosendai.wixsite.com/mysite) in Tohoku University with supports from Daiki Sato and Koki Kido. 
