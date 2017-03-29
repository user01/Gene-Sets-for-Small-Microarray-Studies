
clean:
	-$(RM) results/*.tsv
	-$(RM) results/*.csv
	-find ./plots -name "*.png" -print0 | xargs -0 rm

clean-set-results:
	-$(RM) results/**/*.set.results.*.tsv

test:
	python -m unittest discover
