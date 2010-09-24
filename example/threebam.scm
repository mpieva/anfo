(load 'anfo.scm)
(set-verbosity 'info)

(anfo-run
  (chain
    (concat "foo.anfo" "bar.anfo")
    (filter-score slope: 7.5 intercept: 20))

  (chain
    (only-genomes "hg18")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/hg18.map -o hg18.bam -"))

  (chain
    (only-genomes "pt2")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/pt2.map -o pt2.bam -"))

  (chain
    (only-genomes "hcscca")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/hcscca.map -o hcscca.bam -")))

