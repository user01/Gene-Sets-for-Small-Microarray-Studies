
results/results_all.tsv: results/cluster_em_pca_standard_2.tsv
	Rscript read_results.R

results/cluster_em_pca_standard_2.tsv: results/dimreduced_pca_standard_2.tsv
	Rscript cluster_em.R --name pca_standard_2 --clusters 8 --plot --results

results/dimreduced_pca_standard_2.tsv: results/gene_data_vs_cell_type.tsv
	Rscript dimreduction_pca.R --dimensions 2 --name standard

results/gene_data_vs_cell_type.tsv:
	Rscript load_data.R

clean:
	-$(RM) results/*.tsv
	-$(RM) plots/*.png
