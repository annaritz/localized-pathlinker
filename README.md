# localized-pathlinker
Integrating Protein Localization with Automated Signaling Pathway Reconstruction.

For a signaling pathway of interest, this method identifies paths that connect pathway-specific receptors to pathway-specific transcriptional regulators (TRs) within a weighted protein-protein interactome. We incorporate information about the cellular localization of the interacting proteins to ensure that the individual proteins are localized in cellular compartments involved in signaling transduction and that the protein-protein interactions are spatially coherent with the signaling flow. We model signaling to start at a membrane-bound protein receptor and to be transmitted downstream via a succession of protein-protein interactions within the cytosol to end at a transcription regulator inside the nucleus.

## Installation Instructions
* Localized PathLinker was tested on Python 2.7.9 and requires the following python package(s):
  - NetworkX 1.9.1.
  
* The proteins localization information is needed. Download the "Integrated protein-protein interaction dataset" for the "H. sapiens" from the <a href="http://comppi.linkgroup.hu/downloads">ComPPI databasse</a>.

## Input Files
* NETWORK - A tab-delimited file with one directed interaction per line. Each line should have 3 columns: tail, head, and weight. Edges are directed from tail to head. This represents the interactome. The updated version of the PathLinker interactome, PLNet<sub>2</sub>, can be found in the Data folder under the name "PathLinker_2018_human-ppi-weighted-cap0_75.txt".

* NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first column is the node name and the second is the node type (receptor of TR). You can find an example file for the Alpah6Beta4Integrin pathway in the Data folder under the name "Alpha6Beta4Integrin-nodes.txt"

* COM_PPI - A tab-delimited file with one undirected interaction per line. This represents interactions with localization information. You need to download this file (see Installation Instructions above).

* NODE_LOC_SCORES - A tab-delimited file giving the localization scores per protein for the cellular compartments "ExtMem", "Cytosol", and "Nucleus". You can find an example file derived from the ComPPI database in the Data folder under the name "Protein_Localization_Scores.txt".

## Output Files
* Pathway_paths - A tab-delimited file for the ranked k-shortest paths produced by PathLinker with ties broken. Each path will have two scores: a reconstruction score and a signaling score. Paths with tied reconstruction scores will be re-prioritized by the signaling score.

* Pathway_ranked_edges - A tab-delimiteed file for the edges within the reconstructed paths. Each edge will be given the order of the path within which it appeared for the first time.

## Toy Example
* Assuming your current working directory has the path "Localized_PathLinker", run the following command to apply the "LocPL" technique for the Alpha6Beta4Integrin pathway.

python2.7 /Localized_PathLinker/Loc_PL_run.py -k 20000 -o /Localized_PathLinker/Alpha6Beta4Integrin --write-paths /Localized_PathLinker/2018-03-12-human-ppi-weighted-cap0_75.txt /Localized_PathLinker/Alpha6Beta4Integrin-nodes.txt /Localized_PathLinker/comppi.txt /Localized_PathLinker/Protein_Localization_Scores.txt

## References

This work is currently under review.
The original PathLinker publication is described here: XXX
