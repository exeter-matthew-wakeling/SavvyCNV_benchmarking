# SavvyCNV_benchmarking
Additional files for benchmarking of CNV callers

This repository holds files used for the benchmarking of CNV callers.

## AnalyseCnvs

This is a java program used to calculate the sensitivity and specificity of CNV calling, given a truth set, CNV caller results, and a list of common CNVs. The software can find the optimum filtering parameters for CNVs in the results to produce the maximum true calls and the minimum false calls. The command-line arguments are:
1. -decon - The CNV caller results file is expected to be in the format produced by DeCoN. The file is expected to be a tab-separated file with the sample name in column 2, the chromosome in column 11, the start position in column 9, the end position in column 10, and either "deletion" or "duplication" in column 7. Recommended optimisation arguments are "-optimise 8 -optimise 13 -optimise 0", which correspond to the number of active bins in the CNV, the Bayes Factor, and the Bayes Factor divided by the CNV size.
2. -cnvkit - The CNV caller results file is expected to be in the format produced by CnvKit. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and the log read depth anomaly in column 6. Recommended optimisation arguments are "-optimise 9 -optimise 10 -optimise 0", which correspond to the number of active bins, the quality score, and the quality score divided by the CNV size.
3. -excavator2 - The CNV caller results file is expected to be in the format produced by Excavator2. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and the log read depth anomaly in column 5. Recommended optimisation arguments are "-optimise 9 -optimise 0", which correspond to the CNV probability, and the CNV size.
4. -copywriter - the CNV caller results file is expected to be in the format produced by ConvertCopywriterResults (described below). The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 2, the start position in column 3, the end position in column 4, and either "Deletion" or "Duplication" in column 5. Recommended optimisation arguments are "-optimise 6 -optimise 9 -optimise 10", which correspond to the number of active bins, the absolute value of the log read depth anomaly, and the absolute value of the log read depth anomaly multiplied by the active bin count.
5. -truth - Specifies the location of the truth set. This is a file listing the CNVs that are truly present, as determined by another method, in the same format as the truth set for the ICR96 data set. The file is expected to be a tab-separated file with the sample name in column 1, the chromosome in column 8, the start position in column 9, the end position in column 10, and "Deletion" or "Duplication" in column 6.
6. -optimise (number) - This option can be specified multiple times. It instructs the software to choose a value for the column number specified and filter all the rows in the input file to have a value greater than or equal to this, in such a way as to maximise the specificity and sensitivity. All possible combinations are tried, and the best specificity is printed along with the filter values chosen, for all possible sensitivity values.

Unless otherwise configured, the input file is assumed to be output from SavvyCNV, or from GATK gCNV. The file is expected to be a tab-separated file with the chromosome in column 1, the start location in column 2, the end location in column 3, either "Deletion" or "Duplication" in column 4, and the sample name in column 10. For SavvyCNV, recommended optimisation arguments are "-optimise 7 -optimise 8 -optimise 5", which correspond to the total quality score, the number of active bins, and the quality score divided by the number of active bins. For GATK gCNV, the recommended optimisation arguments are "-optimise 5 -optimise 8 -optimise 6".

The software reads the file "common_cnvs" which should have a list of the locations of common CNVs, which are excluded from the analysis.

If no optimisation options are supplied, the software will output a single row with two tab-separated columns. These are:
1. The number of CNVs that were correctly found.
2. The number of calls in the input file that do not overlap with a true CNV or a common CNV.

If optimisation options are supplied, then multiple rows will be output as many as the number of CNVs that were correctly found. Each row will have the following tab-separated columns:
1. The number of CNVs that were correctly found, after certain filtering values are chosen - this number counts from 1 to the maximum.
2. The lowest possible number of calls in the input file that pass the filters and do not overlap with a true CNV or a common CNV.
3. Multiple columns, giving the values of the chosen filters.

## CalculatePrecisionRecallCurve

Calculates the best precision-recall curve that can be constructed from a set of separate tests. The input is read from the standard input, and is expected to be tab-separated lines. Each line is a separate test. The number of true positive detections and false positives are read from the line. The software takes the following arguments:
1. The first argument is the column number of the true positives.
2. The second argument is the column number of the false positives.
3. The third argument is total number of condition positives (true positives plus false negatives).
4. The fourth argument is optional - if it is not present then the software will emit just the input lines that are on the best precision-recall curve. If the fourth argument is present, then we recommend a value of 1, and the software will interpolate lines between these points. The lines will be curved on a precision-recall graph, but would be straight lines on a ROC curve (but a ROC curve cannot be created without having a well-defined number of true negatives).

The output is a set of tab-separated lines where the first column is the precision and the second column is the recall. Lines that represent an actual input test will have the whole input line appended after these two columns.

## ConvertCopywriterResults

This software converts the results produced by CopywriteR for multiple samples into the result expected by AnalyseCnvs (described above). The software extracts the CNVs from the results, converts chromosomes "23" and "24" into "X" and "Y", and adjusts the expected dosage for male samples. The software takes a single argument, which is the bin size, which is used to find the CopywriteR results. The output is written to the stdout. The CopywriteR output is expected to be in the directory "analysis", containing one directory for each sample, named as the sample name. Each sample directory should have a file called "results_<bin size>" which is processed. The software reads a text file called "male_samples" which is a list of sample names that are male samples, and should have the X chromosome dosage adjusted.

## ProcessGatkGcnvResults

This software converts the results produced by GATK gCNV for multiple samples into the result expected by AnalyseCnvs (described above). GATK gCNV does not emit CNV calls with a start and end position like other CNV callers - rather it produces a list of genomic locations (provided to it) with the probabilities of the different copy numbers for each location. This software converts this into a list of CNVs with a start and end position, by merging together multiple locations that agree on CNV status. The software takes two arguments - the first is the name of a directory containing GATK gCNV results. The second argument is "male" or "female", which adjusts the normal copy number over the X chromosome. The score produced for a CNV is the sum of the probability values over the run of genomic locations, compared to the probability values for the normal copy number. The output is written to the stdout.

## Other files

### Scripts to benchmark tNGS data

| File name | Description |
|------|-------|
| tngs_savvycnv_combinations | This runs SavvyCNV on all the available tNGS data by calling tngs_savvycnv_call_cnvs. The results are placed in files named results_* in the savvycnv directory, which should be created before this script is run. SavvyCNV is run with a bin size of 150000, 200000, 250000, and 300000, with a transition probability of 0.1, 0.03, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, and 0.0000000001, and with a min_reads parameter of 20, 25, 30, 35, 40, 45, and 50. All combinations are run. |
| tngs_savvycnv_call_cnvs | This script is called by tngs_savvycnv_combinations. It runs SavvyCNV six times, to analyse tNGS samples segregated by sex and by sequencing panel. |
| tngs_process_all_savvycnv | This script runs AnalyseCnvs on the output from SavvyCNV. For each file called savvycnv/results_*, it runs AnalyseCnvs three times to find the sensitivity and specificity for all cnvs, large cnvs, and very large cnvs. The results are placed in the savvycnv directory, and named analyse_*, large_analyse_*, and xlarge_analyse_*. |
| tngs_run_all_gatk_gcnv | This runs GATK gCNV on all the available tNGS data by calling tngs_do_gatk_gcnv. GATK gCNV is run with a bin size of 150000, 200000, 250000, and 300000. |
| tngs_do_gatk_gcnv | This script runs GATK gCNV on tNGS data for one particular sex, panel, and bin size. GATK gCNV takes a list of bin locations as an argument, so the list of bins with sufficient coverage was generated by SavvyCNV using the appropriate bin size. As GATK gCNV took so long to run, only the bin list produced by SavvyCNV with a min_reads parameter of 30 were used - this is unlikely to have affected the results significantly. |
| tngs_process_all_gatk_gcnv | This script collects the results from GATK gCNV using ProcessGatkGcnvResults into files named gatk_cnv/results_*, and then analyses them with AnalyseCnvs as with tngs_process_all_savvycnv. |
| tngs_run_all_decon | This runs DeCON on all the available tNGS data by calling tngs_do_decon. DeCON is run with a bin size of 150000, 200000, 250000, and 300000, and with a min_reads parameter of 20, 25, 30, 35, 40, 45, and 50. It runs DeCON six times for each combination to analyse tNGS samples segregated by sex and sequencing panel. |
| tngs_do_decon | This runs DeCON for a particular bin size, min_reads parameter, sex, and sequencing panel. Since DeCON takes a list of bin locations as an argument, the list of bins with sufficient coverage was generated by SavvyCNV using the appropriate bin size and min_reads parameter. The read depth statistics are first analysed and stored, and then the CNVs are called with a transition probability of 0.3, 0.1, 0.03, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, and 0.0000000001. |
| tngs_process_all_decon | This script collects the results from DeCON into files named decon/results_* and then analyses them as with tngs_process_all_savvycnv. |
| tngs_run_all_cnvkit | This script runs CnvKit on all the available tNGS data by calling the six tngs_do_cnvkit_* scripts for the six combinations of sex and sequencing panel. CnvKit is run with bin sizes of 150000, 200000, 250000, 300000, 400000, 500000, 700000, 1000000, 1500000, 2000000, 3000000, and 4000000. |
| tngs_do_cnvkit_* | These scripts run CnvKit for the appropriate samples with the given bin size. CnvKit generates its own bin locations from the bin size, with different size bins for on-target and off-target data. |
| tngs_process_all_cnvkit | This script collects the results from CnvKit into files named cnvkit/results_* and then analyses them as with tngs_process_all_savvycnv. |
| tngs_run_copywriter | This script runs CopywriteR on a single sample, with bin sizes of 50kbp, 100kbp, 200kbp, 500kbp, 1Mbp, and 2Mbp. CopywriteR generates its own bin locations from the bin size, after excluding all data in targeted areas. |
| tngs_convert_copywriter | This script converts the CopywriteR results of multiple samples into a single results_* file in the copywriter directory, using ConvertCopywriterResults. |
| tngs_process_all_copywriter | This script analyses the CNV calls produced by CopywriteR, as with the tngs_process_all_savvycnv script. CopywriteR emits "normal" sections of the genome as well as abnormal copy number CNVs, so this script pre-filters the results to remove sections with an abs(log(dosage)) <= 0.3. |
| tngs_gather_excavator2 | This script gathers the results from Excavator2 into excavator2/results_* files. Excavator2 is run with bin sizes of 150000, 200000, 250000, and 300000. |
| tngs_process_all_excavator2 | This script analyses the CNV calls produced by Excavator2, as with tngs_process_all_savvycnv. |

### Scripts to benchmark exome data

| File name | Description |
|------|-------|
| exome_savvycnv_combinations | This runs SavvyCNV on all the available exome data by calling exome_savvycnv_call_cnvs. The results are placed in files named results_* in the savvycnv directory, which should be created before this script is run. SavvyCNV is run with a bin size of 6000, 7000, 8000, 9000, 10000, 12000, 15000, 20000, 25000, 30000, 40000, and 50000, with a transition probability of 0.1, 0.03, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, and 0.0000000001, and with a min_reads parameter of 20, 25, 30, 35, 40, 45, and 50. All combinations are run. |
| exome_savvycnv_call_cnvs | This script is called by exome_savvycnv_combinations. It runs SavvyCNV eight times, to analyse exome samples segregated by sex and by sequencing panel. |
| exome_process_all_savvycnv | This script runs AnalyseCnvs on the output from SavvyCNV. For each file called savvycnv/results_*, it runs AnalyseCnvs three times to find the sensitivity and specificity for all cnvs and large cnvs. The results are placed in the savvycnv directory, and named analyse_* and large_analyse_*. |
| exome_run_all_gatk_gcnv | This runs GATK gCNV on all the available exome data by calling exome_do_gatk_gcnv. GATK gCNV is run with a bin size of 10000, 12000, 15000, 20000, 25000, 30000, 40000, and 50000. |
| exome_do_gatk_gcnv | This script runs GATK gCNV on exome data for one particular sex, panel, and bin size. GATK gCNV takes a list of bin locations as an argument, so the list of bins with sufficient coverage was generated by SavvyCNV using the appropriate bin size. As GATK gCNV took so long to run, only the bin list produced by SavvyCNV with a min_reads parameter of 30 were used - this is unlikely to have affected the results significantly. |
| exome_process_all_gatk_gcnv | This script collects the results from GATK gCNV using ProcessGatkGcnvResults into files named gatk_cnv/results_*, and then analyses them with AnalyseCnvs as with exome_process_all_savvycnv. |
| exome_run_all_decon | This runs DeCON on all the available exome data by calling exome_do_decon. DeCON is run with a bin size of 10000, 12000, 15000, 20000, 25000, 30000, 40000, and 50000, and with a min_reads parameter of 20, 25, 30, 35, 40, 45, and 50. It runs DeCON eight times for each combination to analyse exome samples segregated by sex and sequencing panel. |
| exome_do_decon | This runs DeCON for a particular bin size, min_reads parameter, sex, and sequencing panel. Since DeCON takes a list of bin locations as an argument, the list of bins with sufficient coverage was generated by SavvyCNV using the appropriate bin size and min_reads parameter. The read depth statistics are first analysed and stored, and then the CNVs are called with a transition probability of 0.3, 0.1, 0.03, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, and 0.0000000001. |
| exome_process_all_decon | This script collects the results from DeCON into files named decon/results_* and then analyses them as with exome_process_all_savvycnv. |
| exome_run_all_cnvkit | This script runs CnvKit on all the available exome data by calling the eight exome_do_cnvkit_* scripts for the eight combinations of sex and sequencing panel. CnvKit is run with bin sizes of 10000, 12000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 100000, 150000, 200000, 300000, and 500000. |
| exome_do_cnvkit_* | These scripts run CnvKit for the appropriate samples with the given bin size. CnvKit generates its own bin locations from the bin size, with different size bins for on-target and off-target data. |
| exome_process_all_cnvkit | This script collects the results from CnvKit into files named cnvkit/results_* and then analyses them as with exome_process_all_savvycnv. |
| exome_run_copywriter | This script runs CopywriteR on a single sample, with bin sizes of 20kbp, 35kbp, 50kbp, 100kbp, 200kbp, 500kbp, 1Mbp, and 2Mbp. CopywriteR generates its own bin locations from the bin size, after excluding all data in targeted areas. |
| exome_convert_copywriter | This script converts the CopywriteR results of multiple samples into a single results_* file in the copywriter directory, using ConvertCopywriterResults. |
| exome_process_all_copywriter | This script analyses the CNV calls produced by CopywriteR, as with the exome_process_all_savvycnv script. CopywriteR emits "normal" sections of the genome as well as abnormal copy number CNVs, so this script pre-filters the results to remove sections with an abs(log(dosage)) <= 0.3. |
| exome_gather_excavator2 | This script gathers the results from Excavator2 into excavator2/results_* files. Excavator2 is run with bin sizes of 10000, 12000, 15000, 20000, 25000, 30000, 40000, and 50000. |
| exome_process_all_excavator2 | This script analyses the CNV calls produced by Excavator2, as with exome_process_all_savvycnv. |
