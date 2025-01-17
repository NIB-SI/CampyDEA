# _Campylobacter jejuni_ biofilm differential expression analysis

This is a repository of Linux commands and R script used for the differential expression analysis of the _C. jenuni_ dataset produced by Manca Volk (ORCID) at Biotechnical Faculty, University of Ljubljana. Dr Marko Petek, NIB ([ORCID](https://orcid.org/0000-0003-3644-7827)) and Manca Volk performed the bioinformatic analysis together. The project was supervised by Dr Špela Baebler, NIB ([ORCID](https://orcid.org/0000-0003-4776-7164)) and Anja Klančnik ([ORCID](https://orcid.org/0000-0003-1632-5785)).

## Content description
- Campy_centrifuge_commands.txt -- taxonomic classification of reads to detect sample contamination
- Campy_STAR_commands.txt -- reference index build and RNA-seq read mapping
- Campy_DEA.R -- differential expression analysis R script

## Bioinformatic software and packages used
- fastqc
- multiqc
- centrifuge
- gffread
- STAR
- samtools
- pavian
- edgeR
- limma
