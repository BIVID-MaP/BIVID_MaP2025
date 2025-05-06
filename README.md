# BIVID_MaP2025

The software for distribution of sequencing reads in SAM file and caluculationg deletion rate for each variant.
![The image of variant-specific interaction profiling](./images/250505_Github_BIVID_MaP.jpg)
# Installation 
```bash
git clone https://github.com/EmiMiyashita-126/BIVID_MaP2025/
```

# Content

# Requirement

・Julia (ver. 1.7.1)

# Usage
## Input

Prepare follwing two arugments: <br>
  `--fasta_path` : Path of the FASTA file containing the reference and mutation sequences<br>
  `--sam_dir` : Path of the directroy including SAM files mapped to reference sequences<br>


## Output

・Divided SAM files for each variant<br>
・Unmapped SAM files containing reads with bases in variant positions deleted<br>
・Base call table containing the number of bases deleted from each base

## File Instructions
・The input FASTA file contains sequence names and their corresponding DNA sequences; the reference sequence is labelled 	`<Seqname>_Ref`.
```text
#./Demo/Input_file/test_FASTA_G4I8.txt
>G4I8_Ref
GAGATGTCTGGCGCAGACATCTCAAATTCAGCGCTTTGGTGGTGGAATGGTGCTATGTGGGCTGAAAAACAAATCGGGCTTCGGTCCGGTTC
>G4I8_Mut
GAGATGTCTGGCGCAGACATCTCAAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAACAAATCGGGCTTCGGTCCGGTTC
```

#Demo
Excuted by MacbookAir (M3)
```bash
cd BIVID_MaP2025
##Output divided SAM files for each variant and base call tables for calculation of deleted reads
julia ./Variant_deletion_profiling.jl --fasta_path ./Demo/Input_file/test_Input_FASTA_G4I8.txt --sam_dir ./Demo/Input_file/input_sam

```
















