#!/bin/bash


alignment_file="test/strap/ZACN.Homo_sapiens.filter2.fasta"
trait_file="test/strap/BodyMass_kg_permulation.cfg"
simulated_traits_file="test/strap/phylorest.tab"
output_file_unfiltered="test/strap/unfiltered.boot"
output_file_filtered="test/strap/filtered.boot"
alignment_format="phylip-relaxed"


python ct bootstrap -a $alignment_file -t $trait_file -s $simulated_traits_file -o $output_file_unfiltered --fmt $alignment_format
python ct bootstrap -a $alignment_file -t $trait_file -s $simulated_traits_file -o $output_file_filtered --fmt $alignment_format --max_miss 0
