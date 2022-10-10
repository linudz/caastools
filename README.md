# caastools v 1.0

## 1. Introduction to CAAStools

Amino acid substitutions that are consistent with phenotypic variation indicate that the gene product is potentially involved in the genetic determination of the trait. We define these cases as Convergent Amino Acid Substitutions (CAAS).

It is possible to retrieve CAAS by scanning Multiple Sequence Alignments (MSA) of orthologous proteins or translated nucleotides. We will isolate those positions in which we can verify that species with diverging trait values converge to different amino acids. A very simple way to do this is to define two groups of species, isolate those positions in which they won’t share any amino acids, and test the statistical significance of this association.

Although the implementation of this strategy can be easily achieved through simple scripting, its scaling to proteome-size requires some additional effort in terms of code optimization. Also, CAAS analysis is usually validated through bootstrap-based approaches. All these operations together are computationally expensive, especially when brought to proteome-wide scale.

In recent years, our group worked on optimizing our in-house scripts for CAAS discovery and validation. CAAStools is the result of these efforts. A small suite of bioinformatics tools that allow the user to identify and validate CAAS analysis on MSA of orthologous proteins.

### 1.1 Suite overview

CAAStools is a collection of 3 bioinformatics tools written in Python 3.9.4. Figure 1 resumes the functioning of each tool. Globally, CAAStools relies on three main pieces of information that are required in different formats (see the input specification section of those paragraphs discussing each single tool). The discovery tool detects CAAS from a single amino acid MSA (protein or translated nucleotide). The resample tools elaborates virtual phenotype groups through different strategies (brownian motion simulation, random sorting of species or phylogeny-restricted simulation). The result of this tool can be submitted to the bootstrap tool that runs a bootstrap CAAS analysis on a single MSA.

![Figure 1]()
