#!/bin/bash


# annotate the genome with AUGUSTUS
./Augustus/src/augustus --species=fly ./Dsim_genome.fasta > ./Output_Annotations/Dsim_annotations.gff
