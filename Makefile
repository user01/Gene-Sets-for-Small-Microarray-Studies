
results/expression_matrix.tsv:
	Rscript load_data.R

clean:
	rm results/*.tsv
