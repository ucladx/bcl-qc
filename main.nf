#!/usr/bin/env nextflow

process align {
  input:
  path fastq_list
  val sample_id
  
  script:
  if( assay == 'pcp' )
    """
    t_coffee -in $sequences > out_file
    """
  else if( assay == 'ces' )
    """
    mafft --anysymbol --parttree --quiet $sequences > out_file
    """
  else 
    error "Invalid alignment mode: ${mode}"
}
  
workflow {
  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola') | sayHello | view
}
