
results/cluster_pca_results.tsv: results/dimreduced_matrix_pca_1.2.tsv
	Rscript cluster_em.R results/dimreduced_matrix_pca_1.2.tsv

results/dimreduced_matrix_pca_1.2.tsv: results/expression_matrix.tsv
	Rscript dimreduction_pca.R

results/expression_matrix.tsv:
	Rscript load_data.R

clean:
	rm results/*.tsv
