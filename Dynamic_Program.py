
from optparse import OptionParser
import sys
from math import log, exp

def Read_Localizatin_Scores(file):

    NODE_LOC_SCORES_file = open(file, 'r')
    Nodes_Loc_Score = {}
    while True:
        line = NODE_LOC_SCORES_file.readline()
        if not line:
            break  # break the while loop. This happens at the end of the file.
        if line == '' or line[0] == '#' or line == '\n':
            continue
        items = line.rstrip().split('\t')
        protein = items[0]
        Nodes_Loc_Score[protein] = {}
        Nodes_Loc_Score[protein]['1'] = max(float(items[1]), 0.01)
        Nodes_Loc_Score[protein]['1'] = -log(Nodes_Loc_Score[protein]['1'])
        Nodes_Loc_Score[protein]['t'] = max(float(items[2]), 0.01)
        Nodes_Loc_Score[protein]['t'] = -log(Nodes_Loc_Score[protein]['t'])
        Nodes_Loc_Score[protein]['T'] = max(float(items[3]), 0.01)
        Nodes_Loc_Score[protein]['T'] = -log(Nodes_Loc_Score[protein]['T'])
    Loc_Keys = Nodes_Loc_Score[protein].keys()  # ['1', 't', 'T']

    return Nodes_Loc_Score, Loc_Keys

def Read_Tied_Paths(file):
    tied_paths_file = open(file, 'r')
    Paths_tied = []
    for line in tied_paths_file:
        if line == '' or line[0] == '#' or line == '\n':
            continue
        items = line.rstrip().split('\t')
        nodes = items[2].rstrip().split('|')
        path_len = len(nodes) # Number of nodes
        if (path_len <= 2): # remove one-edge paths
            continue
        path_temp = [int(items[0])].__add__([float(items[1])]).__add__([nodes])
        Paths_tied.append(path_temp)

    return Paths_tied


def DP_Breaking_Ties(Paths_tied, Nodes_Loc_Score, Loc_Keys):

    Paths_untied = []
    for path in Paths_tied:
        path_nodes = path[2]
        path_len = len(path_nodes)
        Prob_Table = {}
        for i in range(len(Loc_Keys)):
            Prob_Table[Loc_Keys[i]] = [1000] * path_len  # "1000" approximates a probability of zero when log transformed.
        Prob_Table['1'][0] = Nodes_Loc_Score[path_nodes[0]]['1']  # Initialize the first node with its probability to  be in
                                                                  # the ExtMem compartment. The probability of the other
                                                                  # compartments (Cyt, Nuc, and/or Mt) will remain zero.
        for i in range(1, len(path_nodes)):
            Prob_Table['1'][i] = Prob_Table['1'][i-1] + Nodes_Loc_Score[path_nodes[i]]['1']
            
            route1 = Prob_Table['1'][i-1] + Nodes_Loc_Score[path_nodes[i]]['1']
            route2 = Prob_Table['1'][i-1] + Nodes_Loc_Score[path_nodes[i-1]]['t']
            route3 = Prob_Table['t'][i-1]
            Prob_Table['t'][i] = min(route1, route2, route3) + Nodes_Loc_Score[path_nodes[i]]['t']

            route1 = Prob_Table['t'][i-1] + Nodes_Loc_Score[path_nodes[i]]['t']
            route2 = Prob_Table['t'][i-1] + Nodes_Loc_Score[path_nodes[i-1]]['T']
            route3 = Prob_Table['T'][i-1]
            Prob_Table['T'][i] = min(route1, route2, route3) + Nodes_Loc_Score[path_nodes[i]]['T']

        path_score = exp(-Prob_Table['T'][i]) # score of the last protein inside the nucleus.
        path_temp = path.__add__([path_score])
        Paths_untied.append(path_temp)
    Paths_untied = sorted(Paths_untied, key=lambda x: (x[1], x[3]), reverse=True)  # Nested sorting: within paths with equal costs

    return Paths_untied

def print_paths_edges(file, Paths_untied):
    Paths_Str = '%s_untied-paths.txt' % (file)
    Paths_file = open(Paths_Str, 'w')
    Edges_Str = '%s_untied-ranked-edges.txt' % (file)
    Edges_file = open(Edges_Str, 'w')
    Paths_file.write('#New_k\tOld_k\tPath_Length\tPath\tSignaling_Score\n')
    Edges_file.write('#Tail\tHead\tk\n')
    edges = set()
    for k, path in enumerate(Paths_untied, 1):
        nodes = path[2]
        path_len = len(nodes) - 1  # Number of edges
        for i in range(path_len):
            edge = tuple([nodes[i], nodes[i + 1]])
            if edge in edges:
                continue
            edges.add(edge)
            Edges_file.write('%s\t%s\t%d\n' % (nodes[i], nodes[i+1], k))
        Paths_file.write('%d\t%d\t%f\t%s\t%f\n' % (k, path[0], path[1],'|'.join(nodes), path[3]))
    Paths_file.close()
    Edges_file.close()

    print('\nUntied paths are are in: %s' % (Paths_Str))
    print('\nRanked edges are are in: %s' % (Edges_Str))

    return

def main(args):
    usage = '''
Dynamic_Program.py [options] PL_Paths NODE_LOC_SCORES
REQUIRED arguments:
    PL_Paths - A tab-delimited file, where each line is for one path. Each line should
        have three entries. The first one is the path order (k). The second entry is
        the path cost computed by the original PathLinker technique. The third entry
        is the path nodes.

    NODE_LOC_SCORES - A tab-delimited file containing the proteins localizatin socres.
    
'''
    parser = OptionParser(usage=usage)

    parser.add_option('-o', '--output', type='string', default='out_', metavar='STR', \
                      help='A string to prepend to all output files. (default="out")')

    # Parse the command line arguments
    (opts, args) = parser.parse_args()

    # Get the required arguments
    num_req_args = 2
    if len(args) != num_req_args:
        parser.print_help()
        sys.exit('\nERROR: Dynamic_Program.py requires %d positional arguments, %d given.' % (num_req_args, len(args)))
    PL_Paths = args[0]
    localization_scores_file = args[1]

    Nodes_Loc_Score, Loc_Keys = Read_Localizatin_Scores(localization_scores_file)
    Paths_tied = Read_Tied_Paths(PL_Paths)
    Paths_untied = DP_Breaking_Ties(Paths_tied, Nodes_Loc_Score, Loc_Keys)
    print_paths_edges(opts.output, Paths_untied)


if __name__=='__main__':
    main(sys.argv)
