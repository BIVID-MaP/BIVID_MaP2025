# BIVID_MaP2025

The software for distribution of sequencing reads in SAM file and caluculationg deletion rate for each variant.

# Installation 
```bash
git clone https://github.com/…/REPO.git
```

# Content

# Requirement

・Julia (ver. 1.7.1)

# Usage
## Input

・FASTA file including reference sequence and mutant sequences<br>
・SAM file mapped to reference sequence

## Output

・Divided SAM files for each variant<br>
・Unmapped SAM files containing reads with bases in variant positions deleted<br>
・Base call table containing the number of bases deleted from each base

## File Instructions
・The input FASTA file contains sequence names and their corresponding DNA sequences; the reference sequence is labelled 	`<Seqname>_Ref`, while each mutant sequence is labelled 	`<Seqname>_Mut`.












