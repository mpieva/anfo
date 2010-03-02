(require 'anfo)
(require 'pp)
(set-verbosity! 'info)

(pp
  (anfo-run
    "foo.anfo"
    (chain (filter-score slope: 3.5 intercept: 20)
	   (only-genomes '(hg18 pt2))
	   (add-alns)
	   (tee
	     (stats)
	     (mismatches)
	     (divergence primary: 'hg18
			 secondary: 'pt2
			 chop: 4)))))
(newline)

