# _Campylobacter jejuni_ biofilm differential expression analysis (DEA)

This is a repository of Linux commands and R script used for the differential expression analysis of the _C. jenuni_ dataset produced by Manca Volk (ORCID) at Biotechnical Faculty, University of Ljubljana (UNI-LJ BF). Dr Marko Petek, NIB ([ORCID](https://orcid.org/0000-0003-3644-7827)) and Manca Volk performed the bioinformatic analysis together. The project was supervised by Dr Špela Baebler, NIB ([ORCID](https://orcid.org/0000-0003-4776-7164)) and Anja Klančnik, UNI-LJ BF ([ORCID](https://orcid.org/0000-0003-1632-5785)).

## Content description
Run following scripts in consecutive order:
- Campy_centrifuge_commands.txt -- taxonomic classification of reads to detect sample contamination
- Campy_STAR_commands.txt -- reference index build and RNA-seq read mapping
- Campy_DEA.R -- differential expression analysis R script

## Bioinformatic software and packages used
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://github.com/MultiQC/MultiQC)
- [Centrifuge](https://github.com/DaehwanKimLab/centrifuge)
- [GFFRead](https://github.com/gpertea/gffread)
- [STAR](https://github.com/alexdobin/STAR)
- [samtools](https://github.com/samtools/samtools)
- [Pavian](https://github.com/fbreitwieser/pavian)
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)

## R package versions
stringr_1.4.0
edgeR_3.36.0
limma_3.50.3 
