#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:16:23 2020

@author: ligk2e
"""

import os

k = 9   # specify kmer lenght

### Run bowtie_build
bowtie_build_path = "/usr/local/bowtie/1.0.0/bin/bowtie-build"
bowtie_path = "/usr/local/bowtie/1.0.0/bin/bowtie"
#filename = "AE017334.fasta"
#genome_index = "AE017334"
filename = "Bacillus_anthracis_str_sterne.ASM816v1.dna.chromosome.Chromosome.fa"
genome_index = "sterne"

command_build = bowtie_build_path + ' ' + filename + ' ' + genome_index
#os.system(command_build)

### read genome
genome = {}
with open(filename,'r') as f:
    for line in f:
        if line.startswith('>'):
            # do some spliting, here we don't need to do that
            genome[line.rstrip('\n')]=''
            line_use = line.strip('\n')
        else:
            genome[line_use] += line.rstrip('\n')
#print(genome)           
line = '>Chromosome dna:chromosome chromosome:ASM816v1:Chromosome:1:5228663:1 REF' 
genome_seq = genome[line]     
length_genome = len(genome_seq)

### splice genome into kmer
with open('kmer_seq.fasta','w+') as f1:
    
    i = 0
    while i < length_genome - k + 1:
        kmer = genome_seq[i:i+k:1]
        kmer_tandem = kmer * 5
        f1.write('>kmer_3tamdem{}\n'.format(str(i)))
        f1.write(kmer_tandem+'\n')
        i += 1

command_align = bowtie_path + ' ' + genome_index + ' ' +'-f' + ' ' + 'kmer_seq.fasta' + ' ' + '-m 1 -n 0 > sterne_kmer9_5tandem_output'
os.system(command_align)
    
    