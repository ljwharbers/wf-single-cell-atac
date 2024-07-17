#!/bin/bash -euo pipefail
python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
python -c "import sklearn; print(f'scikit-learn,{sklearn.__version__}')" >> versions.txt
fastcat --version | sed 's/^/fastcat,/' >> versions.txt
minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
gffread --version | sed 's/^/gffread,/' >> versions.txt
seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
stringtie --version | sed 's/^/stringtie,/' >> versions.txt
gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
