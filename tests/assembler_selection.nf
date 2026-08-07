#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { assemblersForReadType } from '../lib/workflow_helpers'

workflow {
  ['short', 'nanopore', 'pacbio', 'hifi'].each { read_type ->
    def assemblers = assemblersForReadType(read_type).join(',')
    println "ASSEMBLERS\t${read_type}\t${assemblers}"
  }
}
