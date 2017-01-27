
results/cluster_pca_results.tsv: results/dimreduced_matrix_pca_1.2.tsv
	Rscript cluster_em.R results/dimreduced_matrix_pca_1.2.tsv

results/dimreduced_matrix_pca_1.2.tsv: results/gene_data_vs_cell_type.tsv
	Rscript dimreduction_pca.R

results/gene_data_vs_cell_type.tsv:
	Rscript load_data.R

clean:
	rm results/*.tsv
