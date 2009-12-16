(load 'anfo.scm)
(set-verbosity 'info)

(anfo-run
  (chain
    (concat ; "foo.anfo" "bar.anfo")
      "/mnt/454/NeandertalGenome/PhaseIV/by_library/SL10.anfo"
      "/mnt/454/NeandertalGenome/PhaseIV/by_library/SL11.anfo")
    (filter-score 7.5 20 '()))

  (chain
    (only-genome "hg18")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/hg18.map -o hg18.bam -"))

  (chain
    (only-genome "pt2")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/pt2.map -o pt2.bam -"))

  (chain
    (only-genome "HCSCCA")
    require-hit
    (write-sam "| samtools view -bt /mnt/454/SolexaMPI/PhaseIV/hcscca.map -o hcscca.bam -"))
  )

