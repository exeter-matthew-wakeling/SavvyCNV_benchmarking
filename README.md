# SavvyCNV_benchmarking
Additional files for benchmarking of CNV callers

This repository holds files used for the benchmarking of CNV callers.

## AnalyseCnvs

This is a java program used to calculate the sensitivity and specificity of CNV calling, given a truth set, CNV caller results, and a list of common CNVs. The software can find the optimum filtering parameters for CNVs in the results to produce the maximum true calls and the minimum false calls. The command-line arguments are:
1. -decon - The CNV caller results file is expected to be in the format produced by DeCoN. The file is expected to be a tab-separated file with the sample name in column 2, the chromosome in column 11, the start position in column 9, the end position in column 10, and either "deletion" or "duplication" in column 7. Recommended optimisation arguments are "-optimise 8 -optimise 13 -optimise 0", which correspond to the number of active bins in the CNV, the Bayes Factor, and the Bayes Factor divided by the CNV size.
2. -cnvkit - The CNV caller results file is expected to be in the format produced by CnvKit. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and the log read depth anomaly in column 6. Recommended optimisation arguments are "-optimise 9 -optimise 10 -optimise 0", which correspond to the number of active bins, the quality score, and the quality score divided by the CNV size.
3. -excavator2 - The CNV caller results file is expected to be in the format produced by Excavator2. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and the log read depth anomaly in column 5. Recommended optimisation arguments are "-optimise 9 -optimise 0", which correspond to the CNV probability, and the CNV size.
4. -copywriter - the CNV caller results file is expected to be in the format produced by CopywriteR. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and either "Deletion" or "Duplication" in column 5. Recommended optimisation arguments are "-optimise 6 -optimise 9 -optimise 10", which correspond to the number of active bins, the absolute value of the log read depth anomaly, and the absolute value of the log read depth anomaly multiplied by the active bin count.
5. -truth - Specifies the location of the truth set. This is a file listing the CNVs that are truly present, as determined by another method, in the same format as the truth set for the ICR96 data set. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 8, the start position in column 9, the end position in column 10, and "Deletion" or "Duplication" in column 6.
6. -optimise (number) - This option can be specified multiple times. It instructs the software to choose a value for the column number specified and filter all the rows in the input file to have a value greater than or equal to this, in such a way as to maximise the specificity and sensitivity. All possible combinations are tried, and the best specificity is printed along with the filter values chosen, for all possible sensitivity values.

Unless otherwise configured, the input file is assumed to be output from SavvyCNV, or from GATK gCNV. The file is expected to be a tab-separated file with the chromosome in column 1, the start location in column 2, the end location in column 3, either "Deletion" or "Duplication" in column 4, and the sample name in column 10. For SavvyCNV, recommended optimisation arguments are "-optimise 7 -optimise 8 -optimise 5", which correspond to the total quality score, the number of active bins, and the quality score divided by the number of active bins. For GATK gCNV, the recommended optimisation arguments are "-optimise 5 -optimise 8 -optimise 6".

If no optimisation options are supplied, the software will output a single row with two tab-separated columns. These are:
1. The number of CNVs that were correctly found.
2. The number of calls in the input file that do not overlap with a true CNV or a common CNV.

If optimisation options are supplied, then multiple rows will be output as many as the number of CNVs that were correctly found. Each row will have the following tab-separated columns:
1. The number of CNVs that were correctly found, after certain filtering values are chosen - this number counts from 1 to the maximum.
2. The lowest possible number of calls in the input file that pass the filters and do not overlap with a true CNV or a common CNV.
3. Multiple columns, giving the values of the chosen filters.
