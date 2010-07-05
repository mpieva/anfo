(require 'anfo)
(require 'pp)
(set-verbosity! 'info)

(define (run f)
  (cons f 
        (anfo-run 
          ; (read-file f)
          (apply concat (map read-file (glob f)))
          (chain (tee (stats)
                      (chain (only-genomes '(hg18 pt2))
                             (filter-mapq 60)
                             (filter-total-score slope: 2 intercept: 20 genomes: (list "hg18" "pt2"))
                             (require-hit)
                             (add-alns)
                             (divergence primary: 'hg18 secondary: 'pt2 chop: 4)))))))

(pp (map run (command-line-args)))
(newline)
