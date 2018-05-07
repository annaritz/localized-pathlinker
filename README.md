# localized-pathlinker
Integrating Protein Localization with Automated Signaling Pathway Reconstruction.

For a signaling pathway of interest, this method identifies paths that connect pathway-specific receptors to pathway-specific transcriptional regulators (TRs) within a weighted protein-protein interactome. We incorporate information about the cellular localization of the interacting proteins to ensure that the individual proteins are localized in cellular compartments involved in signaling transduction and that the protein-protein interactions are spatially coherent with the signaling flow. We model signaling to start at a membrane-bound protein receptor and to be transmitted downstream via a succession of protein-protein interactions within the cytosol to end at a transcription regulator inside the nucleus.

## Installation Instructions
* Localized PathLinker was tested on Python 2.7.9 and requires the following python package(s):
  - <a href="https://networkx.github.io/">NetworkX 1.9.1</a>.
  - <a href="https://github.com/Murali-group/PathLinker">PathLinker</a>. Add PathLinker to your environment PATH variable or copy it to your current working directory.

 
## Input Files

These files are positional arguments (the order matters).
* NETWORK - A tab-delimited file that represents the directed, weighted interactome.  Each line contains has to have at least 3 columns: tail, head, and weight. Edges are directed from tail to head. The updated version of the PathLinker interactome, PLNet<sub>2</sub>, can be found in the Data folder under the name "PathLinker_2018_human-ppi-weighted-cap0_75.txt". You need to unzip this file first.

* NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first column is the node name and the second is the node type (receptor, TR, or none). You can find an example file for the Alpah6Beta4Integrin pathway in the Data folder under the name "Alpha6Beta4Integrin-nodes.txt"

* COM_PPI - A tab-delimited file that represents interactions with localization information. We used the predictions from the ComPPI database; you need to download the "Integrated protein-protein interaction dataset" for the "H. sapiens" from the <a href="http://comppi.linkgroup.hu/downloads">ComPPI databasse</a>, and then unzip it. The downloaded file contains one undirected interaction per line. Each line has 12 columns. We use the first (an interaction protein), fifth (the second protein of the interactions), and ninth (interaction score) columns.

* NODE_LOC_SCORES - A tab-delimited file giving the localization scores per protein for the cellular compartments "ExtMem", "Cytosol", and "Nucleus". You can find an example file derived from the ComPPI database in the Data folder under the name "Protein_Localization_Scores.txt". (ARCOMMENT: how did you make this file? Was this file also downloaded from ComPPI?) (ARCOMMENT: after doing this, I would include isntructions about getting hte original file and parsing it.)

## Other Arguments

These arguments may be specified in any order. (ARCOMMENT: note if they are optional or required.)

* -k 

* -o

* --write-paths

(ARCOMMENT: when you run Loc_PL_run.py -h, you should get a help menu. This can also be used here.)

## Output Files
* Pathway_paths - A tab-delimited file for the ranked k-shortest paths produced by PathLinker with ties broken. Each path will have two scores: a reconstruction score and a signaling score. Paths with tied reconstruction scores will be re-prioritized by the signaling score.

* Pathway_ranked_edges - A tab-delimited file for the edges within the reconstructed paths. Each edge will be given the order of the path within which it appeared for the first time.

## Toy Example
This example will compute 20,000 paths for the Alpha6Beta4 Integrin pathway and use the ComPPI information to reweight ties.  Move to the cloned directory (e.g. localized-pathlinker/) and call the following code in python2.7:

```
python Loc_PL_run.py -k 20000 -o example.out --write-paths Data/PathLinker_2018_human-ppi-weighted-cap0_75.txt Data/Alpha6Beta4Integrin-nodes.txt  <COMPPI_INTERACTIONS> Data/Protein_Localization_Scores.txt
```

where <COMPPI_INTERACTIONS> is the interactions file downloaded from ComPPI. 


## References

* This work is currently under review. It is an extension of the method described in this paper:

Ritz A, Poirel CL, Tegge AN, Sharp N, Simmons K, Powell A, Kale SD, and Murali TM, <a href="http://www.nature.com/articles/npjsba20162">Pathways on Demand: Automated Reconstruction of Human Signaling Networks</a>, *npj Systems Biology and Applications*, 2016,2:16002.


* The protein localization information comes from this paper:

Veres DV, Gyurkó DM, Thaler B, Szalay KZ, Fazekas D, Korcsmáros T, and Csermely P, <a href="https://academic.oup.com/nar/article/43/D1/D485/2435307">ComPPI: a cellular compartment-specific database for protein–protein interaction network analysis</a>, *Nucleic Acids Research*, 2015, 43, D1, D485–D493.
