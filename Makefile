
results/results_all.csv: results/cluster_em_pca_standard_2_8c.tsv results/cluster_em_pca_standard_2_12c.tsv results/cluster_em_pca_standard_8_8c.tsv
	Rscript read_results.R

results/cluster_em_pca_standard_2_8c.tsv: results/dimreduced_pca_standard_2.tsv
	Rscript cluster_em.R --name pca_standard_2 --clusters 8 --plot --results

results/cluster_em_pca_standard_2_12c.tsv: results/dimreduced_pca_standard_2.tsv
	Rscript cluster_em.R --name pca_standard_2 --clusters 12 --plot --results

results/cluster_em_pca_standard_8_8c.tsv: results/dimreduced_pca_standard_8.tsv
	Rscript cluster_em.R --name pca_standard_8 --clusters 8 --results


results/dimreduced_pca_standard_2.tsv: results/gene_data_vs_cell_type.tsv
	Rscript dimreduction_pca.R --dimensions 2 --name standard

results/dimreduced_pca_standard_8.tsv: results/gene_data_vs_cell_type.tsv
	Rscript dimreduction_pca.R --dimensions 8 --name standard

results/gene_data_vs_cell_type.tsv:
	Rscript load_data.R

clean:
	-$(RM) results/*.tsv
	-$(RM) plots/*.png
