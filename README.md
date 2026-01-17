# Design and optimization of an integrated computational and machine learning framework to detect functional variants and elucidate their contributions to biological pathways
For the first part of my undergraduate thesis, completed from August 2025 to May 2026 (currently ongoing) at the Komplex Systems Laboratory, Ahmedabad University under the supervision of Professor Krishna Swamy, I am developing a method for the design and optimisation of a computational framework that integrates the analysis of WGS data, including variant calling and Copy Number Variation (CNV) determination. The second part involves the development of an ML ensemble comprising a Naïve Bayes (NB) classifier, a Support Vector Machine (SVM) and an Artificial Neural Network (ANN) that can be used to identify functional variants and predict their possible roles in a biological pathway. The objective is to construct a pipeline that can be reliably trained and tested on prokaryotic as well as eukaryotic data.

WGS reads were aligned to the reference genomes using bwa-mem2. After filtering for duplicate reads and low-quality bases, GATK’s HaplotypeCaller, SelectVariants and VariantFiltration tools were used to identify SNPs and indels. Variant annotation was performed with SnpEff and SnpSift. The Copy Number Variation of the genome was determined with Control-FREEC. 
Each step of the pipeline was run on Stepwell, the university's HPC cluster, which Slurm to schedule and run jobs.

The machine learning pipeline has been built on the code used in [this study](https://www.nature.com/articles/s41467-018-05807-z).
The NB Classifier and SVM will be built in R, which has also been used to preprocess the input data and convert it into a format that is suitable for training ML models. Python will be used to build the ANN.

Datasets from [this paper](https://link.springer.com/article/10.1038/s44320-025-00136-y) were used to test the finished ML pipeline.
