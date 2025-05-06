# BIVID_MaP2025

The software for distribution of sequencing reads in SAM file and caluculationg deletion rate for each variant.

# Installation 
```bash
git clone https://github.com/EmiMiyashita-126/BIVID_MaP2025/
```

# Content

# Requirement

・Julia (ver. 1.7.1)

# Usage
## Input

・FASTA file including reference sequence and mutant sequences<br>
```text
#./Demo/Input_file/test_FASTA_G4I8.txt
>G4I8_Ref
GAGATGTCTGGCGCAGACATCTCAAATTCAGCGCTTTGGTGGTGGAATGGTGCTATGTGGGCTGAAAAACAAATCGGGCTTCGGTCCGGTTC
>G4I8_Mut
GAGATGTCTGGCGCAGACATCTCAAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAACAAATCGGGCTTCGGTCCGGTTC
```
・SAM file mapped to reference sequence
```text
#./Demo/Input_file/test_TGIRT.sam
@SQ	SN:G4I8_Ref	LN:92
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem /Users/emimiyashita/Desktop/事務書類/2025etal/analysis/Figure_SI_Run252_G4I8_G4timedependency/Fastq/reference/I8_ref/I8_ref /Users/emimiyashita/Desktop/事務書類/2025etal/analysis/Figure_SI_Run252_G4I8_G4timedependency/Fastq/Fastp_output/S19_S19_L001_R1_001_fastp.fastq /Users/emimiyashita/Desktop/事務書類/2025etal/analysis/Figure_SI_Run252_G4I8_G4timedependency/Fastq/Fastp_output/S19_S19_L001_R2_001_fastp.fastq
@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.12	CL:samtools view -h -q 20 /Users/emimiyashita/Desktop/事務書類/2025etal/analysis/Figure_SI_Run252_G4I8_G4timedependency/Fastq/SAM/S19_S19_removed.sam
FS10003367:8:BTR67813-1623:1:1101:1090:1000	97	G4I8_Ref	24	60	32M	=	24	32	AAATTCAGCGCTTTGGTGGTGGAATGGTGCTA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:32	MC:Z:32M	AS:i:32	XS:i:0
FS10003367:8:BTR67813-1623:1:1101:1090:1000	145	G4I8_Ref	24	60	32M	=	24	-32	AAATTCAGCGCTTTGGTGGTGGAATGGTGCTA	::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:32	MC:Z:32M	AS:i:32	XS:i:0
FS10003367:8:BTR67813-1623:1:1101:1230:1000	99	G4I8_Ref	24	60	46M	=	24	46	AAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	MD:Z:26G19	MC:Z:46M	AS:i:41	XS:i:0
FS10003367:8:BTR67813-1623:1:1101:1230:1000	147	G4I8_Ref	24	60	46M	=	24	-46	AAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAA	FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	MD:Z:26G19	MC:Z:46M	AS:i:41	XS:i:0
FS10003367:8:BTR67813-1623:1:1101:1560:1000	99	G4I8_Ref	24	60	46M	=	24	46	AAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	MD:Z:26G19	MC:Z:46M	AS:i:41	XS:i:0
FS10003367:8:BTR67813-1623:1:1101:1560:1000	147	G4I8_Ref	24	60	46M	=	24	-46	AAATTCAGCGCTTTGGTGGTGGAATGATGCTATGTGGGCTGAAAAA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:1	MD:Z:26G19	MC:Z:46M	AS:i:41	XS:i:0
```

## Output

・Divided SAM files for each variant<br>
・Unmapped SAM files containing reads with bases in variant positions deleted<br>
・Base call table containing the number of bases deleted from each base

## File Instructions
・The input FASTA file contains sequence names and their corresponding DNA sequences; the reference sequence is labelled 	`<Seqname>_Ref`, while each mutant sequence is labelled 	`<Seqname>_Mut`.












