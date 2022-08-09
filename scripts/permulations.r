#!/usr/bin/env Rscript

#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/


# A Convergent Amino Acid Substitution identification 
# and analysis toolbox
# 
# Author:         Fabio Barteri (fabio.barteri@upf.edu)
# 
# Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
#                 Xavier Farré (xfarrer@igtp.cat),
#                 David de Juan (david.juan@upf.edu).
# 
# 
# SCRIPT NAME: permulations.r
# DESCRIPTION: Permulation script from RERconverge
# DEPENDENCIES: modules in modules/simulation folder
# CALLED BY: simulatte.py



# Import RERconverge (with required package install)

if (!require("RERconverge", character.only=T, quietly=T)) {
  require(devtools)
  install_github("nclark-lab/RERconverge", ref="master") 
}

library(RERconverge)
library(tibble)
library(readr)
library(ape)

# Inputs 

args = commandArgs(trailingOnly=TRUE)

tree <- args[1]
config.file <- args[2]
number.of.cycles <- args[3]
selection.strategy <- args[4]
phenotypes <- args[5]
outfile <- args[6]

#tree <- "inputs/traits.speciestree.roooted.210907.nh"
#config.file <-"inputs/BodyMass_kg_permulation.cfg"
#number.of.cycles <- "1000"
#phenotypes <- "inputs/BodyMass_kg_permulation.txt"
#outfile <- "prova.txt"

# Read the tree object.
#imported.tree <- read_file(tree)

tree.o <- read.tree(tree)
trait <- tree.o$tip.label
l <- length(trait)

# Read the config file 

cfg <- read.table(config.file, sep ="\t", header = F)

foreground.df <- subset(cfg, cfg$V2 == "1")
background.df <- subset(cfg, cfg$V2 == "0")

foreground.size <- length(foreground.df$V1)
background.size <- length(background.df$V1)

foreground.species <- foreground.df$V1
background.species <- background.df$V1

# Read the phenotype file

phenotype.df <- read.table(phenotypes, sep = "\t")
foreground.values <- subset(phenotype.df, phenotype.df$V1 %in% foreground.species)$V2
background.values <- subset(phenotype.df, phenotype.df$V1 %in% background.species)$V2


# SEED AND PRUNE
starting.values <- phenotype.df$V2
all.species <- phenotype.df$V1

pruned.tree.o <- drop.tip(tree.o,tree.o$tip.name[-match(all.species, tree.o$tip.name)])

names(starting.values) <- all.species

### SIMULATE
counter = 0

simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))

for (j in 1:as.integer(number.of.cycles)){
  counter = counter + 1
  cycle.tag = paste("b", as.character(counter), sep = "_")
  permulated_phenotype <- simpermvec(starting.values, pruned.tree.o)
  x <- enframe(permulated_phenotype)

  # Select potential foreground and background species
  potential.fg.df <- subset(x, value %in% foreground.values)
  potential.bg.df <- subset(x, value %in% background.values)
  
  if (selection.strategy == "random") {
  
    potential.fg <- potential.fg.df$name
    potential.bg <- potential.bg.df$name

    fg.species <- sample(potential.fg, foreground.size)
    bg.species <- sample(potential.bg, background.size)

  }

  else if (selection.strategy == "inner") {
  
    fg.species <- head(potential.fg.df[order(potential.fg.df$value, decreasing=TRUE), ], foreground.size)$name
    bg.species <- head(potential.bg.df[order(potential.bg.df$value, decreasing=FALSE), ], background.size)$name

  }

  else if (selection.strategy == "edges") {

    fg.species <- head(potential.fg.df[order(potential.fg.df$value, decreasing=FALSE), ], foreground.size)$name
    bg.species <- head(potential.bg.df[order(potential.bg.df$value, decreasing=TRUE), ], background.size)$name

  }

  else {

    paste("Wrong species selection options for permulations. Please select between 'random', 'inner' and 'edges'.")
  
  }

  fg.species.tag = paste(fg.species, collapse=",")
  bg.species.tag = paste(bg.species, collapse=",")
  
  logline = paste("Processing permulation", as.character(counter), "of", as.character(number.of.cycles))
  outline = c(cycle.tag, fg.species.tag, bg.species.tag)
  #simulated.traits.df <- rbind(simulated.traits.df, outline)
  simulated.traits.df[nrow(simulated.traits.df) +1,] <- outline
  
  write(logline, stdout())
}

# Export the result
write.table(simulated.traits.df, sep="\t",col.names=FALSE, row.names = FALSE, file=outfile, quote=FALSE)

