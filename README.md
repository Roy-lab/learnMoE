# learnMoE : GMM clustering
original : SR's original code, pulled from programs/clustering/ on Jan 25 2018

printPredictions : DC's branch to print out predicted expression, per-module results correlations, and overall train/test correlation results (when in CV mode). Based on original. Jan 25 2018
	New files are fold*/pred_test.txt and cc_test.txt.
	run_example.sh provides guidance for running an example.

learnMoE_tab : edited from printPredictions. This one can read in a tab delimited file in the same format as our expression matrices: first column is gene name, subsequent columns are samples. Genes are treated as variables.
