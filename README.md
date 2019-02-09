# localized-pathlinker
Integrating Protein Localization with Automated Signaling Pathway Reconstruction.

For a signaling pathway of interest, this method identifies paths that connect pathway-specific receptors to pathway-specific transcriptional regulators (TRs) within a weighted protein-protein interactome. We incorporate information about the cellular localization of the interacting proteins to ensure that the individual proteins are localized in cellular compartments involved in signaling transduction and that the protein-protein interactions are spatially coherent with the signaling flow. We model signaling to start at a membrane-bound protein receptor and to be transmitted downstream via a succession of protein-protein interactions within the cytosol to end at a transcription regulator inside the nucleus.

## Installation Instructions
* Localized PathLinker was tested on Python 2.7.9 and requires the following python package(s):
  - <a href="https://networkx.github.io/">NetworkX 1.9.1</a>.
  - <a href="https://github.com/Murali-group/PathLinker">PathLinker</a>. Add PathLinker to your environment PATH variable or copy it to your current working directory.

 
## Input Files

These files are positional arguments (the order matters).
* NETWORK - A tab-delimited file that represents the directed, weighted interactome.  Each line has to have at least 3 columns: tail, head, and weight. Edges are directed from tail to head. The updated version of the PathLinker interactome, *PLNet<sub>2</sub>*, can be found in the Data folder under the name "PathLinker_2018_human-ppi-weighted-cap0_75.txt". You need to unzip this file first.

* NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first column is the node name and the second is the node type (receptor, TR, or none). You can find an example file for the Alpah6Beta4Integrin pathway in the Data folder under the name "Alpha6Beta4Integrin-nodes.txt"

* COM_PPI - A tab-delimited file that represents interactions with localization information. We used the predictions from the ComPPI database; you need to download the "Integrated protein-protein interaction dataset" for the "H. sapiens" from the <a href="http://comppi.linkgroup.hu/downloads">ComPPI databasse</a>, and then unzip it. The downloaded file contains one undirected interaction per line. Each line has 12 columns. We use the first (an interaction protein), fifth (the second protein of the interaction), and ninth (interaction score) columns.

* NODE_LOC_SCORES - A tab-delimited file giving the proteins localization scores- one protein per line. Each line has four columns: protein ID, "ExtMem" localization score, "Cytosol" localization score, and "Nucleus" localization score. You can find an example file derived from the ComPPI database in the Data folder under the name "Protein_Localization_Scores.txt". This file was obtained using the python script "Nodes_Scores.py" available with the source codes. Read the required arguments in "Nodes_Scores.py" to know the needed files to derive the "Protein_Localization_Scores.txt" and how to obtain them.

## Other Arguments

These optional arguments may be specified in any order.

* -k: number of computed paths. Default is 100 paths.

* -o: name to proceed output files. Default is 'out'.

* --write-paths: If given, the computed paths are saved in a file, in addition to the ranked edges.


## Output Files
* Pathway_paths - A tab-delimited file for the ranked k-shortest paths produced by PathLinker with ties broken. Each path will have two scores: a reconstruction score and a signaling score. Paths with tied reconstruction scores will be re-prioritized by the signaling score.

* Pathway_ranked_edges - A tab-delimited file for the edges within the reconstructed paths. Each edge will be given the order of the path within which it appeared for the first time.

## Toy Example
This example will compute 20,000 paths for the Alpha6Beta4 Integrin pathway and use the ComPPI information to reweight ties.  Move to the cloned directory (e.g. localized-pathlinker/) and call the following code in python2.7:

```
python Loc_PL_run.py -k 20000 -o example.out --write-paths Data/PathLinker_2018_human-ppi-weighted-cap0_75.txt Data/Alpha6Beta4Integrin-nodes.txt  <COMPPI_INTERACTIONS> Data/Protein_Localization_Scores.txt
```

where <COMPPI_INTERACTIONS> is the interactions file downloaded from ComPPI. 


<!--- ## *PLNet<sub>2</sub>* Interactome -->

We built *PLNet<sub>2</sub>* from both physical molecular interaction data (BioGrid, DIP, InnateDB, IntAct, MINT, PhosphositePlus) and annotated signaling pathway databases (KEGG, NetPath, and SPIKE) [3-7]. *PLNet<sub>2</sub>* contains 17,168 nodes, 40,016 directed regulatory interactions, and 286,250 bidirected physical interactions, totaling 612,516 directed edges. We assigned interaction direction based on evidence of a directed enzymatic reaction (e.g., phosphorylation, dephosphorylation, ubiquitination) from any of the source databases.  Each interaction is supported by one or more types of experimental evidence (e.g. yeast two hybrid or co-immunoprecipitation) that are available as evidence codes from the data sources, and/or the name of the pathway database it is from. Edges are weighted using an evidence-based Bayesian approach that assigns higher confidence to an experiment type/pathway database if it identifies interacting proteins that participate in the same biological process [8]. Given a set *P* of positive edges and a set *N* of negative edges, the method estimates, for each evidence type *t*, a probability that *t* supports positive interactions. These probabilities are then combined for each interaction supported by (potentially multiple) evidence types to produce a final weight. We chose the GO term "regulation of signal transduction" to build a set of positive interactions that are likely related to signaling; this term includes "signal transduction" as a child GO term.
Positives are edges whose nodes are both annotated with this term, and negatives are randomly selected edges whose nodes are not co-annotated to the term. We chose *|N|* = 10 x *|P|* negative edges. To lessen the influence of very highly-weighted edges, we apply a ceiling of 0.75 to all weights [8]. 

## References

* **This work is currently under review. It is an extension of the method described in this paper:**

[1] Ritz A, Poirel CL, Tegge AN, Sharp N, Simmons K, Powell A, Kale SD, and Murali TM, <a href="http://www.nature.com/articles/npjsba20162">Pathways on Demand: Automated Reconstruction of Human Signaling Networks</a>, *npj Systems Biology and Applications*, 2016,2:16002.


* **The protein localization information comes from this paper:**

[2] Veres DV, Gyurkó DM, Thaler B, Szalay KZ, Fazekas D, Korcsmáros T, and Csermely P, <a href="https://academic.oup.com/nar/article/43/D1/D485/2435307">ComPPI: a cellular compartment-specific database for protein–protein interaction network analysis</a>, *Nucleic Acids Research*, 2015, 43, D1, D485–D493.


* **Resources used to build *PLNet<sub>2</sub>*:**

[3] Aranda B *et al.*, <a href="https://www.nature.com/articles/nmeth.1637">PSICQUIC and PSISCORE: accessing and scoring molecular interactions</a>, *Nature Methods*, 2011, 8, 528–529.

[4] Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, and Sullivan M, <a href="https://academic.oup.com/nar/article/40/D1/D261/2903142">PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse</a>, *Nucleic Acids Research*, 2012, 40, D1, D261–D270.

[5] Kandasamy K *et al.*, <a href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-1-r3">NetPath: a public resource of curated signal transduction pathways</a>, *Genome Biology*, 2010, 11, R3.

[6] Kanehisa M, Furumichi M, Tanabe M, Sato Y, and Morishima K, <a href="https://academic.oup.com/nar/article/45/D1/D353/2605697">KEGG: new perspectives on genomes, pathways, diseases and drugs</a>, *Nucleic Acids Research*, 2017, 45, D1, D353–D361.

[7] Paz A *et al.*, <a href="https://academic.oup.com/nar/article/39/suppl_1/D793/2507440">SPIKE: a database of highly curated human signaling pathways</a>, *Nucleic Acids Research*, 2011, 39, Issue suppl_1, D793–D799.

[8] Yeger-Lotem E *et al.*, <a href="https://www.nature.com/articles/ng.337">Bridging high-throughput genetic and transcriptional data reveals cellular responses to alpha-synuclein toxicity</a>, *Nature Genetics*, 2009, 41, 316–323.
