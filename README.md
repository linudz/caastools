
# CAAStools 1.0 - Software documentation.
# 1. Introduction to CAAStools

Amino acid substitutions that are consistent with phenotypic variation indicate that the gene product is potentially involved in the genetic determination of the trait. We define these cases as _Convergent Amino Acid Substitutions_ (CAAS).

It is possible to retrieve CAAS by scanning Multiple Sequence Alignments (MSA) of orthologous proteins or translated nucleotides. We will isolate those positions in which we can verify that species with diverging trait values converge to different amino acids. A very simple way to do this is to define two groups of species, isolate those positions in which they won’t share any amino acids, and test the statistical significance of this association.

Although the implementation of this strategy can be easily achieved through simple scripting, its scaling to proteome-size requires some additional effort in terms of code optimization. Also, CAAS analysis is usually validated through bootstrap-based approaches. All these operations together are computationally expensive, especially when brought to proteome-wide scale.

In recent years, our group worked on optimizing our in-house scripts for CAAS discovery and validation. CAAStools is the result of these efforts. A small suite of bioinformatics tools that allow the user to identify and validate CAAS analysis on MSA of orthologous proteins.


## 1.1 Suite overview 

CAAStools is a collection of 3 bioinformatics tools written in Python 3.9.4. **Figure 1** resumes the functioning of each tool. Globally, CAAStools relies on three main pieces of information that are required in different formats (see the input specification section of those paragraphs discussing each single tool). The **_discovery_** tool detects CAAS from a single amino acid MSA (protein or translated nucleotide). The **_resample_** tools elaborates virtual phenotype groups through different strategies (brownian motion simulation, random sorting of species or phylogeny-restricted resampling). The result of this tool can be submitted to the **_bootstrap_** tool that runs a bootstrap CAAS analysis on a single MSA.

[Link to **Figure 1**](https://figshare.com/articles/figure/CAAStools_Figure/21324306)

The output of the discovery tool consists in a table reporting a list of CAAS associations. A corresponding p-value is calculated as the probability of randomly finding a CAAS association in that position (see our preprint for a full explanation about the statistical testing of CAAS).


# 2 Download, installation and updates.


## 2.1 Availability.

CAAStools is available as a GitHub repository

[https://github.com/linudz/caastools](https://github.com/linudz/caastools).  

The repository can be cloned through command line

`git clone github.com/linudz/caastools`

Also, you can clone it through the GitHub GUI client or download the code as a zipped file through your browser.


## 2.2 Installation.


The file _ct_ in the CAAStools Git folder is an executable python script that will launch the different tools of the suite. To run it system-wide, you can decide to either download the code to the usr/local/bin folder or to add the caastools folder to the $PATH variable.


### Updates.

Please, follow the github.com/linudz/caastools repository for code updates.


## 2.3 Dependencies and software requirements.

CAAStools is written in **Pyhton 3.9.4** and is compatible with any Python 3+ environment. Each tool has its own set of dependencies.


### Discovery and Bootstrap tools.

Biopython 1.79

Scipy (Biopython dependency)

Numpy (Biopython dependency)


### Resample tool

Dendropy 4+

Also, the _Brownian Motion_ _mode_ (`--mode bm`) for trait simulation relies on the simpermvec() function from R library RERconverge ([https://github.com/nclark-lab/RERconverge](https://github.com/nclark-lab/RERconverge)). This library will be requested only by the _resample_ tool when run in _Brownian Motion mode_ (`--mode bm`).

**Warning.** Apple M1 and M2 users might experience issues during RERconverge installation ([https://github.com/nclark-lab/RERconverge/issues/60](https://github.com/nclark-lab/RERconverge/issues/60)).


# 3. Discovery tool


## 3.1 Algorithm overview

CAAStools discovery identifies CAAS from an amino acid Multiple Sequence Alignment file (MSA). It is possible to view the inputs and the options through the `-h --help` command

`ct discovery -h/--help`

CAAStools discovery returns those positions in which amino acids differ between two groups of species. We call these groups **discovery groups**, distinguishing them into foreground (**FG**) and background (**BG**). The program will detect a CAAS in those MSA positions that will meet two conditions. First, the foreground and the background won’t share any amino acids. Second, all the species in at least one group need to share (or _converge_ to) the same amino acid. The two conditions define a set of 4 possible mutation **patterns**. _Table 3.1_ reports a set of possible mutation patterns, indicating which ones are accepted as CAAS, that are enumerated from 1 to 4.


<table>
  <tr>
   <td><strong>FG</strong>
   </td>
   <td><strong>BG</strong>
   </td>
   <td><strong>Difference between groups</strong>
   </td>
   <td><strong>Convergence within…</strong>
   </td>
   <td><strong>Pattern</strong>
   </td>
   <td><strong>Is it CAAS?</strong>
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>W
   </td>
   <td>YES
   </td>
   <td>Both groups
   </td>
   <td>Pattern 1
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>WE
   </td>
   <td>YES
   </td>
   <td>FG
   </td>
   <td>Pattern 2
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>W
   </td>
   <td>YES
   </td>
   <td>BG
   </td>
   <td>Pattern 3
   </td>
   <td>YES
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>WF
   </td>
   <td>YES
   </td>
   <td>None
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>KE
   </td>
   <td>NO
   </td>
   <td>FG
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
  <tr>
   <td>KE
   </td>
   <td>K
   </td>
   <td>NO
   </td>
   <td>BG
   </td>
   <td>-
   </td>
   <td>NO
   </td>
  </tr>
</table>


**Table 3.1** - Mutation patterns.

Note that the **pattern 4** can be included as CAAS by user specification. Through the `--patterns` option, the user will be able to select the number of patterns to include in the output. By default, ct discovery returns the patterns 1,2 and 3.

The user can filter the result based on the maximum number of indels (or gaps, “-”) accepted per position or the maximum species missing in the alignment.

For each CAAS prediction, the program will calculate an empiric p-value that corresponds to the probability of finding the same set of mutational patterns with random species. This probability is calculated through the hypergeometric distribution.

For further details on the CAAStools discovery algorithm, please refer to our preprint on BioRXIV.


## 3.2 Formatting the inputs


### 3.2.1 The configuration file

`-t /--traitfile $config_file`

The **configuration file **or **config** is the file that we’ll use to tell the program which species we are comparing and how they are arranged into the FG and BG. It consists of a simple tab file containing the name of the species and a label indicating the corresponding group (**FG** = 1, **BG **= 0).

A config file is present in the examples/ folder (examples/conifig.tab)

Aotus_griseimembra	0

Avahi_peyrierasi	0

Callibella_humilis	0

Gorilla_beringei	1

Gorilla_gorilla	1

Macaca_thibetana	1

Mandrillus_leucophaeus	1

This config file will instruct the program to create two groups.


<table>
  <tr>
   <td>FG
   </td>
   <td>Gorilla_beringei, Gorilla_gorilla, Macaca_thibetana, Mandrillus_leucophaeus
   </td>
  </tr>
  <tr>
   <td>BG
   </td>
   <td>Aotus_griseimembra, Avahi_peyrierasi, Callibella_humilis
   </td>
  </tr>
</table>


Note that:



* The values are tab-separated and no further space is admitted
* The order of the species is irrelevant.


### 3.2.2 The amino acid MSA

The second fundamental input is the amino acid MSA file. CAAStools relies on Biopython 1.7 to import sequence files. Hence, the accepted formats are the ones specified in Biopython docs ( [https://biopython.org/wiki/AlignIO](https://biopython.org/wiki/AlignIO) ):

**clustal, emboss, fasta, fasta-m10, ig, maf, mauve, msf, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm**

By default, ct discovery will read the MSA as a clustal file. To specify a different format, we’ll need to specify it through the `--fmt` option (e.g. `--fmt phylip-relaxed`). In the examples folder, the examples/MSA directory contains an MSA in different formats.


### 3.2.3 A note on name consistency

CAAStools will associate the sequence in the MSA to the species by name identity. The program will save the name of each species in FG and BG groups in string variables, and will select those sequences in the MSA whose ID field will coincide with one of the species in the config file. **Please, format your config file and MSA in order to match the sequence IDs with the name of the species in the config file**.


## 3.3 Gaps and missing species filtering 

CAAS results can be filtered for a maximum of gaps or missing species. In this, we define as a “gap” the presence of an indel which is indicated with the “-” character. A missing species will be a species that is mentioned in the config file but it is not found in the MSA. This situation can occur when we iterate the analysis over different MSAs with variable coverage.

By default, CAAStools discovery accepts a maximum of n-1 gaps or missing species per group, where n is the size of the group. This means that the program will need the presence of at least one species per group to verify the conditions for CAAS assignment. The user can decide to limit the number of gaps and missing species per group, or to skip those positions in which gaps represent more than a maximum percentage of total symbols (default=50%). The following tables report the different options for gaps and missing species.


### 3.3.1 Filtering for gaps


<table>
  <tr>
   <td>Max background gaps
   </td>
   <td>--max_bg_gaps
   </td>
   <td>Filter by number of gaps in the background.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max foreground gaps
   </td>
   <td>--max_fg_gaps
   </td>
   <td>Filter by number of gaps in the foreground.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max overall gaps
   </td>
   <td>--max_fg_gaps
   </td>
   <td>Filter by total number of gaps
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max gaps per position
   </td>
   <td>--max_gaps_per_position
   </td>
   <td>Filter by number of gaps per position.
   </td>
   <td>Default: 0.5 (50%)
   </td>
  </tr>
</table>



### 3.3.2 Filtering for missing species


<table>
  <tr>
   <td>Max background missing species
   </td>
   <td>--max_bg_miss
   </td>
   <td>Filter by number of missing species in the background.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max foreground missing species
   </td>
   <td>--max_fg_miss
   </td>
   <td>Filter by number of missing species in the foreground.
   </td>
   <td>Default: No filter
   </td>
  </tr>
  <tr>
   <td>Max overall missing species
   </td>
   <td>--max_fg_miss
   </td>
   <td>Filter by total number of missing species
   </td>
   <td>Default: No filter
   </td>
  </tr>
</table>

## 3.4 The output

CAAStools discovery will output a tab file with all the CAAS found in one single MSA.


<table>
  <tr>
   <td><strong>Column</strong>
   </td>
   <td><strong>Header</strong>
   </td>
   <td><strong>Description</strong>
   </td>
  </tr>
  <tr>
   <td>1
   </td>
   <td>Gene
   </td>
   <td>The name of the gene (from MSA filename)
   </td>
  </tr>
  <tr>
   <td>2
   </td>
   <td>Trait
   </td>
   <td>The name of the trait (from binary config file)
   </td>
  </tr>
  <tr>
   <td>3
   </td>
   <td>Position
   </td>
   <td>The position in the MSA (0-based)
   </td>
  </tr>
  <tr>
   <td>4
   </td>
   <td>Substitution
   </td>
   <td>The substitution FG/BG
   </td>
  </tr>
  <tr>
   <td>5
   </td>
   <td>Pvalue
   </td>
   <td>The p-value from hypergeometric distribution
   </td>
  </tr>
  <tr>
   <td>6
   </td>
   <td>Scenario
   </td>
   <td>The mutational pattern (see “patterns” table in<em> 3.1 - Algorithm overview</em>)
   </td>
  </tr>
  <tr>
   <td>7
   </td>
   <td>FFGN
   </td>
   <td>Species <strong>F</strong>ound in <strong>FG</strong>: <strong>N</strong>umber. Number of species found in the FG group (it excludes those ones having an indel). 
   </td>
  </tr>
  <tr>
   <td>8
   </td>
   <td>FBGN
   </td>
   <td>Species <strong>F</strong>ound in <strong>BG</strong>: <strong>N</strong>umber. Number of species found in the BG group (it excludes those ones having an indel). 
   </td>
  </tr>
  <tr>
   <td>9
   </td>
   <td>GFG
   </td>
   <td>Number of <strong>G</strong>aps in the <strong>FG</strong>.
   </td>
  </tr>
  <tr>
   <td>10
   </td>
   <td>GBG
   </td>
   <td>Number of <strong>G</strong>aps in the <strong>BG</strong>.
   </td>
  </tr>
  <tr>
   <td>11
   </td>
   <td>MFG
   </td>
   <td>Number of <strong>M</strong>issing species in the <strong>FG</strong>.
   </td>
  </tr>
  <tr>
   <td>12
   </td>
   <td>MBG
   </td>
   <td>Number of <strong>M</strong>issing species in the <strong>BG</strong>.
   </td>
  </tr>
  <tr>
   <td>13
   </td>
   <td>FFG
   </td>
   <td>List of species <strong>F</strong>ound in the <strong>FG</strong>. Comma-separated.
   </td>
  </tr>
  <tr>
   <td>14
   </td>
   <td>FBG
   </td>
   <td>List of species <strong>F</strong>ound in the <strong>FG</strong>. Comma-separated.
   </td>
  </tr>
  <tr>
   <td>15
   </td>
   <td>MS
   </td>
   <td>List of missing species. Comma-separated.
   </td>
  </tr>
</table>



## 


## 3.5 Examples

Run ct discovery with example alignment (`phylip-relaxed` format, that has to be specified).

`ct discovery -a examples/MSA/primates.msa.pr -t examples/config.tab -o examples/discovery.output.usr.example --fmt phylip-relaxed`


# 4 Resample tool


## 4.1 Algorithm Overview

CAAStools resample (ct resample) elaborates a set of resampled discovery groups for bootstrap analysis. The global options for the tool can be fetched through the help command:

`ct resample -h/--help`

The simulation of discovery groups is propaedeutic to bootstrap analysis

can consist in a simple randomization, a randomization that is restricted to some parts of the phylogeny, or be based on brownian-motion trait evolution simulation. The user will indicate one of these three strategies, the size of the resampled FG and BG groups and the number of simulation cycles. The program outputs a tab file in 


### 4.1.1 Random simulation strategy

`ct resample --mode random`

In this case, the simulation will consist in the bare random sorting of species into a pair of FG/BG discovery groups. 


### 4.1.2 Phylogeny-restricted random simulation strategy.

`ct resample --mode random --limit_by_group $groupfile`

In this case, the simulation is based on the random choice of species, but is limited to the families that are present in a config file provided as a template. A further file, the species file, specifies the composition of the families. The random scooping takes into account the number of groups (or families) present in the template groups and will replicate that composition. For instance, if our template FG group consists of 3 species from group A and 2 species from groupB, the randomisation will follow this pattern. In each cycle, the program scoops 3 random species from group A and 2 random species from group B.


### 4.1.3 Brownian motion based simulation strategy.

`ct resample --mode bm`

This strategy resamples neutral evolution by brownian motion simulation. The FG/BG group size is defined by a template config file. The R function

`simpermvec()`

from R library RERconverge ([https://github.com/nclark-lab/RERconverge](https://github.com/nclark-lab/RERconverge)) will perform a brownian motion simulation of neutral evolution. FG and BG will be defined as the n-th species with higher values and the m-th species with lower values, where n and m are the size of FG and BG respectively.


## 4.2 Inputs per simulation strategy

Each simulation strategy will require a specific set of input files. The following table reports all the inputs that are needed for each strategy.


<table>
  <tr>
   <td><strong>Input File</strong>
   </td>
   <td><strong>Random</strong>
   </td>
   <td><strong>Random (Phylogeny restricted)</strong>
   </td>
   <td><strong>Brownian Motion</strong>
   </td>
  </tr>
  <tr>
   <td>Phylogenetic tree in newick format.
<p>
-p/--phylogeny
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
  </tr>
  <tr>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Config File as a template
<p>
--bytemp
   </td>
   <td>Can be replaced by -f/--fg_size and -b/--bg_size for FG/BG size definition
   </td>
   <td>Mandatory
   </td>
   <td>Mandatory
   </td>
  </tr>
  <tr>
   <td>Trait values
<p>
--traitvalues
   </td>
   <td>NO
   </td>
   <td>NO
   </td>
   <td>Mandatory. It is used by the program to shuffle the
   </td>
  </tr>
</table>

The user will need to provide a phylogenetic tree in Newick format ([https://evolution.genetics.washington.edu/phylip/newicktree.html](https://evolution.genetics.washington.edu/phylip/newicktree.html)) to tell the program on which species base its simulation. The binary configuration file is required for the phylogeny restricted and Brownian motion strategies. The phylogeny restricted strategy limits trait randomisation to some specific families.

The resample tool outputs 1000 resampled traits by default. The user can decide the number of cycles through the --cycles option:

`--cycles 10`

`--cycles 100`

`--cycles 1000`


## 4.3 The output

The output consists in a tab file with three columns.

**Column 1**: The name of the cycle, indicated in the b_numberofcycle format

**Column 2**: The FG species (comma-separated)

**Column 3**: The BG species (comma-separated)

Here’s an example of the first ten lines of a resampled traits output file.

  b_1	Cercopithecus_mitis,Alouatta_palliata,Pan_troglodytes,Colobus_polykomos	Saimiri_sciureus,Alouatta_puruensis,Trachypithecus_crepusculus,Cercocebus_torquatus

  b_2	Macaca_fuscata,Papio_cynocephalus,Cercopithecus_petaurista,Cheracebus_lucifer	Trachypithecus_phayrei,Tarsius_wallacei,Lophocebus_aterrimus,Macaca_silenus

  b_3	Saimiri_oerstedii,Nomascus_concolor,Lemur_catta,Saguinus_oedipus	Pongo_abelii,Indri_indri,Eulemur_albifrons,Eulemur_fulvus

  b_4	Saimiri_ustus,Eulemur_rubriventer,Leontocebus_nigricollis,Macaca_mulatta	Semnopithecus_hypoleucos,Mus_musculus,Otolemur_garnettii,Eulemur_macaco

  b_5	Cacajao_hosomi,Alouatta_puruensis,Saimiri_macrodon,Pithecia_mittermeieri	Ateles_belzebuth,Macaca_maura,Prolemur_simus,Trachypithecus_laotum

  b_6	Eulemur_coronatus,Trachypithecus_geei,Hapalemur_griseus,Prolemur_simus	Eulemur_flavifrons,Macaca_fuscata,Trachypithecus_pileatus,Cercopithecus_mona

  b_7	Cheracebus_lugens,Cercocebus_lunulatus,Hapalemur_occidentalis,Ateles_chamek	Tarsius_lariang,Cercopithecus_neglectus,Cercopithecus_diana,Macaca_thibetana

  b_8	Perodicticus_potto,Pan_paniscus,Cacajao_hosomi,Lepilemur_ankaranensis	Cercocebus_chrysogaster,Papio_anubis,Callimico_goeldii,Plecturocebus_miltoni

  b_9	Macaca_leonina,Cercopithecus_diana,Propithecus_diadema,Macaca_siberu	Galagoides_demidovii,Alouatta_seniculus,Propithecus_coquereli,Plecturocebus_dubius

  b_10	Leontopithecus_rosalia,Papio_cynocephalus,Hylobates_agilis,Mus_musculus	Cheracebus_regulus,Cercopithecus_mitis,Lophocebus_aterrimus,Ateles_marginatus


## 4.4 Examples {#4-4-examples}


### Ex.1 Resampling based on random selection of species {#ex-1-resampling-based-on-random-selection-of-species}

**With input fg/bg size (-f and -b options)**

 `ct resample -p examples/phylogeny.nw -f 5 -b 4 -m random --cycles 500 -o test/resample/random.resampling.tab`

**By template (binary config)**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --cycles 500 -o test/resample/random.resampling.bytemplate.tab`

**Phylogeny restricted (must go by template)**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --limit_by_group test/sp2fam.210727.tab --cycles 500 -o test/resample/random.resampling.bytemplate.tab`


### Ex.2 resampling based on BM 

**Template and trait values mandatory**

`ct resample -p examples/phylogeny.nw --bytemp examples/config.tab -m random --cycles 500  --traitvalues examples/traitvalues.tab -o test/resample/BM.resampling.tab`




# 5 Bootstrap Tool


## 5.1 Algorithm overview 

The bootstrap tool is designed to repeat the CAAS discovery on a large number of discovery groups. In this case, the discovery groups are defined through the output file of the resample tool, in which each line represents a single cycle (see previous paragraph). The program will scan an MSA and will count the number of cycles that return a CAAS in that position.


## 5.2 The inputs

`-s $resampled_trait (output of ct resample)`

`-a $MSA`


## 5.3 Output

The output consists in a tabbed file with three columns

Column 1: Position

Column 2: Number of resamples returning a CAAS in the position

Column 3: Number of cycles

Column 4: Bootstrap value

Column 5: Cycles with positive CAAS


## 5.4 Examples


### Bootstrap from random resampled traits.

`ct bootstrap -s test/resample/random.resampling.tab -a examples/MSA/primates.msa.pr -o examples/random.bootstrap.tab --fmt phylip-relaxed -t test/BodyMass_kg_permulation.cf`

5. License

This software is licensed under GNU General Public License. The kind of license is to be decided with UPF.


# 6. How to cite

The application note for CAAStools is [available as a preprint in bioRxiv](https://doi.org/10.1101/2022.12.14.520422).


# 7. How to contact the development team

For any inquire, please contact Fabio Barteri at Pompeu Fabra University / BBRC [fabio.barteri@upf.edu](mailto:fabio.barteri@upf.edu)

