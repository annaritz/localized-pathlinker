import PathLinker as PL
import ksp_Astar as KSP
import Dynamic_Program as DP

import sys
from optparse import OptionParser, OptionGroup
import networkx as nx


def main(args):
    usage = '''
Localized_PathLinker.py [options] NETWORK NODE_TYPES COM_PPI NODE_LOC_SCORES
REQUIRED arguments:
    NETWORK - A tab-delimited file with one directed interaction per
        line. Each line should have at least 2 columns: tail, head. Edges
        are directed from tail to head. This file can have a third column
        specifying the edge weight, which is required unless the --PageRank
        option is used (see --PageRank help for a note on these weights).
        To run PathLinker on an unweighted graph, set all edge weights
        to 1 in the input network.

    NODE_TYPES - A tab-delimited file denoting nodes as receptors or TRs. The first
        column is the node name, the second is the node type, either 'source'
        (or 'receptor') or 'target' (or 'tr' or 'tf'). Nodes which are neither receptors nor TRs may
        be omitted from this file or may be given a type which is neither 'source'
        nor 'target'.

    COM_PPI - A tab-delimited file with one undirected interaction per line. This file is a cellular
        compartment-sepcific interactome.
    
    NODE_LOC_SCORES - 


'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('-o', '--output', type='string', default='out_', metavar='STR', \
                      help='A string to prepend to all output files. (default="out")')

    parser.add_option('', '--write-paths', action='store_true', default=False, \
                      help='If given, also output a list of paths found by KSP in addition to the ranked edges.')

    parser.add_option('', '--no-log-transform', action='store_true', default=False, \
                      help='Normally input edge weights are log-transformed. This option disables that step.')

    parser.add_option('', '--largest-connected-component', action='store_true', default=False, \
                      help='Run PathLinker on only the largest weakly connected component of the graph. May provide performance speedup.')

    parser.add_option('', '--edge-penalty', type='float', default=1.0, \
                      help='Factor by which to divide every edge weight. The effect of this option is to penalize the score of every path by a factor equal to (the number of edges in the path)^(this factor). (default=1.0)')

    parser.add_option('', '--PL_Weights', type='int', default=1, metavar='INT',
                     help='If (-1), use the weights of the ComPPI interactome.')

    # # Random Walk Group
    # group = OptionGroup(parser, 'Random Walk Options')
    #
    # group.add_option('', '--PageRank', action='store_true', default=False, \
    #                  help='Run the PageRank algorithm to generate edge visitation flux values, which are then used as weights for KSP. A weight column in the network file is not needed if this option is given, as the PageRank visitation fluxes are used for edge weights in KSP. If a weight column is given, these weights are interpreted as a weighted PageRank graph.')
    #
    # group.add_option('-q', '--q-param', action='store', type='float', default=0.5, \
    #                  help='The value of q indicates the probability that the random walker teleports back to a source node during the random walk process. (default=0.5)')
    #
    # group.add_option('-e', '--epsilon', action='store', type='float', default=0.0001, \
    #                  help='A small value used to test for convergence of the iterative implementation of PageRank. (default=0.0001)')
    #
    # group.add_option('', '--max-iters', action='store', type='int', default=500, \
    #                  help='Maximum number of iterations to run the PageRank algorithm. (default=500)')
    #
    # parser.add_option_group(group)

    # k shortest paths Group
    group = OptionGroup(parser, 'k Shortest Paths Options')

    group.add_option('-k', '--k-param', type='int', default=100, \
                     help='The number of shortest paths to find. (default=100)')

    group.add_option('', '--allow-mult-targets', action='store_true', default=False, \
                     help='By default, PathLinker will remove outgoing edges from targets to ensure that there is only one target on each path.  If --allow-mult-targets is specified, these edges are not removed.')

    group.add_option('', '--allow-mult-sources', action='store_true', default=False, \
                     help='By default, PathLinker will remove incoming edges to sources to ensure that there is only one source on each path.  If --allow-mult-sources is specified, these edges are not removed.')

    parser.add_option_group(group)

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    # get the required arguments
    num_req_args = 4
    if len(args) != num_req_args:
        parser.print_help()
        sys.exit('\nERROR: PathLinker.py requires %d positional arguments, %d given.' % (num_req_args, len(args)))

    NETWORK_FILE = args[0]
    NODE_VALUES_FILE = args[1]
    ComPPI_FILE = args[2]
    NODE_LOC_SCORES_FILE = args[3]

    # Read the ComPPI interactome edges
    ComPPI_file = open(ComPPI_FILE, 'r')
    ComPPI_edges = {}
    ComPPI_edges_zero = set()
    ComPPI_file.readline()  # To read the header line
    for line in ComPPI_file:
        if line == '' or line[0] == '#' or line == '\n':
            continue
        items = line.rstrip().split('\t')
        edge_temp = tuple(sorted([items[0], items[4]]))
        ComPPI_edges[edge_temp] = float(items[8])
        if float(items[8]) == 0:
            ComPPI_edges_zero.add(edge_temp)

    ComPPI_edges_Final = set(ComPPI_edges.keys()).difference(ComPPI_edges_zero)
    ComPPI_nodes_Final = set([u for u, v in ComPPI_edges_Final]).union(set([v for u, v in ComPPI_edges_Final]))


    net = None

    # Read the network file
    print('\nReading the network from %s' % (NETWORK_FILE))
    with open(NETWORK_FILE, 'r') as network_file:
        net = PL.readNetworkFile(network_file, False)

    print(nx.info(net))
    # Remove edges that do not have cellular localization information.
    for u,v in net.edges():
        edge_temp = tuple(sorted([u,v]))
        if edge_temp not in ComPPI_edges_Final:
            net.remove_edge(u,v)

    # Print info about the network
    print(nx.info(net))

    # Read the sources and targets on which to run PageRank and KSP
    sources = set()
    targets = set()

    # Read the receptors and TRs file
    print("Reading sources and targets from " + NODE_VALUES_FILE)
    for line in open(NODE_VALUES_FILE, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split('\t')]

        # Skip empty lines and lines beginning with '#' comments
        if line == '':
            continue
        if line[0] == '#':
            continue

        if items[1] in ['source', 'receptor'] and items[0] in ComPPI_nodes_Final:
            sources.add(items[0])
        elif items[1] in ['target', 'tr', 'tf'] and items[0] in ComPPI_nodes_Final:
            targets.add(items[0])

    print('\nRead %d sources and %d targets' % (len(sources), len(targets)))

    # Remove sources and targets that don't appear in the network, and do some sanity checks on sets
    sources = set([s for s in sources if s in net])
    targets = set([t for t in targets if t in net])
    print('\tAfter removing sources and targets that are not in the network: %d sources and %d targets.' % (
        len(sources), len(targets)))
    if len(sources) == 0:
        sys.exit('ERROR: No sources are in the network.')
    if len(targets) == 0:
        sys.exit('ERROR: No targets are in the network.')
    if len(sources.intersection(targets)) > 0:
        sys.exit('ERROR: %d proteins are listed as both a source and target.' % (len(sources.intersection(targets))))

    ## Prepare the network to run KSP

    # Remove improper edges from the sources and targets. This portion
    # must be performed before the log transformation, so that the
    # renormalization within accounts for the probability lost to the
    # removed edges.  These transformations are executed by default;
    # to prevent them, use the opts.allow_mult_sources or opts.allow_mult_targets
    # arguments.
    if not opts.allow_mult_sources:
        PL.modifyGraphForKSP_removeEdgesToSources(net, sources)
    if not opts.allow_mult_targets:
        PL.modifyGraphForKSP_removeEdgesFromTargets(net, targets)

    # Apply the user specified edge penalty
    PL.applyEdgePenalty(net, opts.edge_penalty)

    # Transform the edge weights with a log transformation
    if (not opts.no_log_transform):
        PL.logTransformEdgeWeights(net)

    # Add a super source and super sink. Performed after the
    # transformations so that the edges can be given an additive
    # weight of 0 and thus not affect the resulting path cost.
    PL.modifyGraphForKSP_addSuperSourceSink(net, sources, targets, weightForArtificialEdges=0)

    ## Run the pathfinding algorithm
    print('\nComputing the k=%d shortest simple paths.' % (opts.k_param))
    paths = KSP.k_shortest_paths_yen(net, 'source', 'sink', opts.k_param, weight='ksp_weight')

    if len(paths) == 0:
        sys.exit('\tERROR: Targets are not reachable from the sources.')

    ## Use the results of KSP to rank edges

    # Un-does the logarithmic transformation on the path lengths to
    # make the path length in terms of the original edge weights
    if not opts.no_log_transform:
        paths = PL.undoLogTransformPathLengths(paths)


    # Preparing the paths for the Dynamic Program
    k = 0
    paths_PL = [] # Original paths by PathLinker
    for path in paths:
        pathNodes = [n for n, w in path]
        pathNodes = pathNodes[1:-1]
        path_len = len(pathNodes) # Number of nodes
        if (path_len <= 2): # remove one-edge paths
            continue
        k += 1
        path_cost = path[-1][1]
        path_temp = [k].__add__([path_cost]).__add__([pathNodes])
        paths_PL.append(path_temp)

    Nodes_Loc_Score, Loc_Keys = DP.Read_Localizatin_Scores(NODE_LOC_SCORES_FILE)
    Paths_untied = DP.DP_Breaking_Ties(paths_PL, Nodes_Loc_Score, Loc_Keys)
    DP.print_paths_edges(opts.output, Paths_untied)

    print('\nFinished!')

if __name__ == '__main__':
    main(sys.argv)


