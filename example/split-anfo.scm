(load 'anfo.scm)
(set-verbosity 'info)

(define (seq a b) (if (eq? a b) '() (cons a (seq (+ a 1) b))))

(define genome "hg18")
(define input-file "Vi-hg18.anfo")
(define output-file "Vi-hg18-chr")
(define sequences 
  (map (lambda (c) (string-append "chr" c)) 
       (append (map number->string (seq 1 22)) '("X" "Y" "M"))))

(define (output s)
  (chain 
    (require-hit genome s)
    (write-native (string-append output-file s ".anfo") 99)))

(apply anfo-run input-file (map output sequences))
