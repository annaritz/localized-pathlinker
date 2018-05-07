import sys
from optparse import OptionParser, OptionGroup


def main(args):
    usage = '''
Localized_PathLinker.py NETWORK COM_PPI NODE_LOC_SCORES NODE_LOC_SCORES_out
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per
        line. Each line should have at least 2 columns: tail, head. Edges
        are directed from tail to head. This file can have a third column
        specifying the edge weight, which is required unless the --PageRank
        option is used (see --PageRank help for a note on these weights).
        To run PathLinker on an unweighted graph, set all edge weights
        to 1 in the input network.

    COM_PPI - A tab-delimited file with one undirected interaction per line. This file is a cellular
        compartment-sepcific interactome.

    NODE_LOC_SCORES - A tab-delimited file downloaded from the ComPPI database. Visit the following link:
        "http://comppi.linkgroup.hu/downloads", and choose "Integrated subcellular localization dataset" for 
        "H. sapiens" with "All Localizations". Unzip the file and use it here.

    NODE_LOC_SCORES_out - The output file path and name. The output file will be a tab-delimited file with one node per
        line. Each line has four columns:
            1- Node ID.
            2- Maximum localization score of Extracelluar and Membrane.
            3- Localization score at the Cytosol.
            4- Localization score at the Nucleus.
'''

    parser = OptionParser(usage=usage)
    (opts, args) = parser.parse_args()

    # get the required arguments
    num_req_args = 4
    if len(args) != num_req_args:
        sys.exit('\nERROR: PathLinker.py requires %d positional arguments, %d given.' % (num_req_args, len(args)))

    NETWORK_FILE = args[0]
    ComPPI_FILE = args[1]
    NODE_LOC_SCORES_FILE = args[2]
    NODE_LOC_SCORES_FILE_out = args[3]

    # Read the PLNet2 interactome edges.
    PL_file = open(NETWORK_FILE, 'r')
    PL_edges = set()
    PL_file.readline()  # To read the header line
    for line in PL_file:
        if line == '' or line[0] == '#' or line == '\n':
            continue
        items = line.rstrip().split('\t')
        edge_temp = tuple(sorted([items[0], items[1]]))
        PL_edges.add(edge_temp)

    # Read the ComPPI interactome edges
    ComPPI_file = open(ComPPI_FILE, 'r')
    ComPPI_edges = set()
    ComPPI_edges_zero = set()
    ComPPI_file.readline()  # To read the header line
    for line in ComPPI_file:
        if line == '' or line[0] == '#' or line == '\n':
            continue
        items = line.rstrip().split('\t')
        edge_temp = tuple(sorted([items[0], items[4]]))
        ComPPI_edges.add(edge_temp)
        if float(items[8]) == 0:
            ComPPI_edges_zero.add(edge_temp)

    ComPPI_edges_Final = ComPPI_edges.difference(ComPPI_edges_zero)
    common_edges_Final = ComPPI_edges_Final.intersection(PL_edges)
    common_nodes_Final = set([u for u, v in common_edges_Final]).union(set([v for u, v in common_edges_Final]))

    # Output file
    outfile = open(NODE_LOC_SCORES_FILE_out, 'w')
    outfile.write('#Node\tExtMem\tCytosol\tNucleus\n')

    Loc_file = open(NODE_LOC_SCORES_FILE, 'r')
    Loc_file.readline()
    for line in Loc_file:
        items = line.rstrip().split('\t')
        if items[0] in common_nodes_Final:
            items_temp = items[3]
            items_temp = items_temp.rstrip().split('|')

            Prob_th = 0  # If you need to include localization scores above a certain value.
            CompProb_temp1 = 0
            CompProb_tempt = 0
            CompProb_tempT = 0
            for u in items_temp:
                u_temp = u.rstrip().split(':')
                if float(u_temp[1]) >= Prob_th:
                    if u_temp[0] == 'extracellular' or u_temp[0] == 'membrane':
                        if CompProb_temp1 == 0:  # To get the maximum probability of either 'E' or 'M'
                            CompProb_temp1 = float(u_temp[1])
                        elif CompProb_temp1 < float(u_temp[1]):
                            CompProb_temp1 = float(u_temp[1])
                    elif u_temp[0] == 'cytosol':
                        CompProb_tempt = float(u_temp[1])
                    elif u_temp[0] == 'nucleus':
                        CompProb_tempT = float(u_temp[1])
            outfile.write('%s\t%0.6f\t%0.6f\t%0.6f\n' % (items[0], CompProb_temp1, CompProb_tempt, CompProb_tempT))


if __name__ == '__main__':
    main(sys.argv)