(require 'pp)
(require 'anfo)
(set-verbosity! 'info)

(define maps "/nfs/eva/udo/anfo")
(define oname (car (command-line-args)))
(define patterns (cdr (command-line-args)))

(define (bam-out genome)
  (chain
    (only-genomes genome)
    (require-hit)
    (write-sam (string-append "| samtools view -bt " maps "/" genome
                              ".map -o " oname "-" genome ".bam -"))))

(define (div-filtered slope)
  (chain
    (filter-total-score slope: slope)
    (fmap (lambda (r)
            (cons (cons 'slope slope) r))
          (divergence primary: 'hg18 secondary: 'pt2 chop: 4))))

(define run-it
  (anfo-run 
    (apply concat (map read-file (apply append (map glob patterns))))
    (filter-score)
    (sanitize)
    (sort-pos '("hg18" "pt2") mem: 2000 handles: 512)
    (rmdup)
    (apply tee
           (stats)
           (write-native (string-append oname ".anfo"))
           (chain (add-alns)
                  (apply tee
                         (map div-filtered '(7.5 3.5 2 1))))
           (map bam-out '("hg18" "pt2" "hcscca")))))

(with-output-to-file
  (string-append oname ".txt")
  (lambda () (pp run-it) (newline)))

