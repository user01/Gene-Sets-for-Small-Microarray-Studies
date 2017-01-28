
results/cluster_pca_results.tsv: results/dimreduced_pca_standard_2.tsv
	Rscript cluster_em.R results/dimreduced_pca_standard_2.tsv

results/dimreduced_pca_standard_2.tsv: results/gene_data_vs_cell_type.tsv
	Rscript dimreduction_pca.R --dimensions 2 --name standard

results/gene_data_vs_cell_type.tsv:
	Rscript load_data.R

clean:
	-$(RM) results/*.tsv
	-$(RM) plots/*.png
