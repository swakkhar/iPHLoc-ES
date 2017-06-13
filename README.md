# iPHLoc-ES

iPHLoc-ES: Identification of Bacteriophage Protein Locations using Evolutionary and Structural Features

## installation:

no installation needed. following prerequisite are required to run the software. just download the files in a folder or download the iPHLoc-ES.tar.gz file and extract it into your local directory 

1. install python 3.4 or latest version
2. install following python packages:
	- numpy
	- pandas
	- pickle
	- sklearn

## run the software:

Use the following commands from the terminal:

python3.4 predict.py PSSM SPD3 OPTION


Usage:/home/swakkhar/Dropbox/Protein-subcelluar localization/iPHLoc-ES/predict.py PSSM SPD3 OPTION

PSSM: PSSM file generated from PSI-BLAST

SPD3: SPD3 file generated from SPIDER3

OPTION = 0 for identification of cellular phages

OPTION = 1 for prediction of phage location

