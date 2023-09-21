#!/bin/bash

output_file_unfiltered="test/strap/unfiltered.boot"
output_file_filtered="test/strap/filtered.boot"


# Step 1. Select a species to be modified

species="Tarsius_wallacei"              # Now this is Tarxius_wallacei

# Step 2. Select the cycles where this species appears          LIST 1

cat test/strap/phylorest.tab | grep $species | awk '{print $1}' | sort | uniq > test/strap/twal.cycles 

# Step 5. Fetch all the cycles from unfiltered bootstrap        LIST 2

cat $output_file_unfiltered | 
awk -F "\t" '{print $4}' |
sort | uniq |
tr "," "\n" |
sort | uniq > test/strap/unfiltered.cycles

# Step 6. Fetch all the cycles from filtered bootstrap          LIST 3

cat $output_file_filtered | 
awk -F "\t" '{print $4}' |
sort | uniq |
tr "," "\n" |
sort | uniq > test/strap/filtered.cycles
