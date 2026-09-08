.PHONY: verify verify-core verify-symbolic paper clean

verify: verify-core verify-symbolic

verify-core:
	python verification/verify_phylogenetic_Z17.py
	python verification/audit_Z17_independent.py

verify-symbolic:
	python verification/audit_Z17_symbolic.py

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error phylogenetic_Z17_counterexample_EN.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error phylogenetic_Z17_counterexample_EN.tex

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc
