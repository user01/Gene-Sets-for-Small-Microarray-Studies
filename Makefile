
clean:
	-$(RM) results/*.tsv
	-$(RM) pairwise/*.tsv
	-$(RM) results/*.csv
	-find ./plots -name "*.png" -print0 | xargs -0 rm
