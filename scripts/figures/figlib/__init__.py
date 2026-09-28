"""
Shared code for the manuscript figures (FIGURE_2.py ... FIGURE_7.py).

  paths       where pipeline outputs live (params.local.yml, results/ layout)
  cohort      groups, display labels, colours and markers (the Figure 1 palette)
  style       matplotlib settings, panel labels, saving, tables and reports
  annotation  imprinting centre, BP clusters, segmental duplications, genes
  stats       rank tests, permutation tests, bootstrap, BH (numpy only)
  readers     VCF, HiFiCNV, pbsv and FASTA-index readers
  modbam      per-molecule CpG calls from MM/ML-tagged BAMs (needs pysam)
"""
