import pandas as pd
import numpy
import networkx

def getConnectionToTF(model, allTFs, affectedNodes):
    g = networkx.from_pandas_edgelist(model, 'source', 'target', create_using=networkx.DiGraph())
    allNodes = numpy.array(list(g.nodes))
    includedTFs = numpy.intersect1d(allTFs, allNodes)

    connectedToTF = numpy.isin(affectedNodes, includedTFs)
    for i in range(len(affectedNodes)):
        if affectedNodes[i] in allNodes:
            for tf in includedTFs:
                if networkx.algorithms.shortest_paths.generic.has_path(g, affectedNodes[i], tf):
                    connectedToTF[i] = True
                    break

    return connectedToTF

def getConnectionToLigand(model, allLigands, affectedNodes):
    g = networkx.from_pandas_edgelist(model, 'target', 'source', create_using=networkx.DiGraph())
    allNodes = numpy.array(list(g.nodes))
    includedLigands = numpy.intersect1d(allLigands, allNodes)
    connectedToLigand = numpy.isin(affectedNodes, allLigands)
    for i in range(len(affectedNodes)):
        if affectedNodes[i] in allNodes:
            for ligand in includedLigands:
                if networkx.algorithms.shortest_paths.generic.has_path(g, affectedNodes[i], ligand):
                    connectedToLigand[i] = True
                    break
    return connectedToLigand

#def trimDeadEnds(model, allTFs, allLigands):
def trimDeadEnds(model, allTFs,maintain):
    allNodes = numpy.union1d(model['source'].values, model['target'].values)
    allMaintain = numpy.union1d(maintain['source'].values,maintain['target'].values)
    maintainedNodes = numpy.isin(allNodes , allMaintain)

    connectedToTF = getConnectionToTF(model, allTFs, allNodes)
    #connectedToLigand = getConnectionToLigand(model, allLigands, allNodes)
    #connectedToBoth = numpy.logical_and(connectedToTF, connectedToLigand)
    #disconectedNodes = allNodes[connectedToBoth == False]
    disconectedNodes = allNodes[numpy.logical_and(connectedToTF == False , maintainedNodes==False)]

    dissconectedSources = numpy.isin(model.source.values, disconectedNodes)
    disconnectedTargets = numpy.isin(model.target.values, disconectedNodes)
    disconnectedEdges = numpy.logical_or(dissconectedSources, disconnectedTargets)
    model = model.loc[disconnectedEdges==False,:]
    return model

#def trimSelfConnecting(model, allTFs, allLigands):
def trimSelfConnecting(model, allTFs,maintain):
    allMaintain = numpy.union1d(maintain['source'].values, maintain['target'].values)
    lastSize = numpy.inf
    curentSize = model.shape[0]
    while curentSize<lastSize:
        lastSize = curentSize
        sources, counts = numpy.unique(model['source'], return_counts=True)
        onlyOneInput = sources[counts == 1]
        targets, counts = numpy.unique(model['target'], return_counts=True)
        onlyOneOutput = targets[counts == 1]
        overlap = numpy.intersect1d(onlyOneInput, onlyOneOutput)
        #Exclude ligands and TFs
        #overlap = numpy.setdiff1d(overlap, allLigands)
        overlap = numpy.setdiff1d(overlap, allTFs)
        selfLoop = numpy.full(len(overlap), False, dtype=bool)
        for i in range(len(selfLoop)):
            curSource = model.loc[numpy.isin(model['target'], overlap[i]), 'source'].values[0]
            curTarget = model.loc[numpy.isin(model['source'], overlap[i]), 'target'].values[0]
            selfLoop[i] = curSource == curTarget

        affectedProteins = overlap[selfLoop]

        affectedInteractions = numpy.logical_or(numpy.isin(model['source'], affectedProteins), numpy.isin(model['target'], affectedProteins))
        maintainedInteractions= numpy.logical_or(numpy.isin(model['source'], allMaintain), numpy.isin(model['target'], allMaintain))
        affectedInteractions = numpy.logical_and(numpy.logical_not(maintainedInteractions),affectedInteractions)
        model = model.loc[affectedInteractions==False,:]
        curentSize = model.shape[0]

    return model

def subsetOnSource(df, coreSources):
    dfFilter = numpy.full(df.shape[0], False, dtype=bool)
    for i in range(len(dfFilter)):
        dfFilter[i] = len(numpy.intersect1d(df.iloc[i,:]['source'].split(';'), coreSources))>0
    df = df.loc[dfFilter,:].copy()
    return df

def MaxSubgraph(model):
    allNodes = numpy.union1d(model['source'].values, model['target'].values)
    g = networkx.from_pandas_edgelist(model, 'target', 'source', create_using=networkx.DiGraph())
    ug = g.to_undirected()
    sub_graphs = [ug.subgraph(c).copy() for c in networkx.connected_components(ug)]
    max_sg=networkx.Graph()
    n=0
    for i,sg in enumerate(sub_graphs):
        if n<=sg.number_of_nodes():
            n = sg.number_of_nodes()
            max_sg = sg.copy()

    conectedSources = numpy.isin(model.source.values, numpy.array(list(max_sg.nodes())))
    connectedTargets = numpy.isin(model.target.values, numpy.array(list(max_sg.nodes())))
    connectedEdges = numpy.logical_or(conectedSources, connectedTargets)
    model = model.loc[connectedEdges == True, :]
    return model

coreSources = ['KEGG', 'InnateDB','SIGNOR']


#Load and trim PKN
PKN = pd.read_csv('../preprocessing/preprocessed_data/FilteredOmnipath_v2.tsv', sep='\t', low_memory=False, index_col=0)
PKNFull = PKN.copy()

TFgene = pd.read_csv('../results/filtered_shrnas_tf_activities.tsv', sep='\t', low_memory=False, index_col=0)
MaintainInteractions = pd.read_csv('../preprocessing/preprocessed_data/MaintainInter.tsv', sep='\t', low_memory=False, index_col=0)
targetdTFs = pd.read_csv('../preprocessing/preprocessed_data/targetd_tfs.tsv', sep='\t', low_memory=False, index_col=0)
allTFs = numpy.intersect1d(TFgene.columns.values, PKN['target'].values)
allTFs = numpy.union1d(allTFs,targetdTFs['Entry'].values)

PKN = trimDeadEnds(PKN, allTFs,MaintainInteractions)
PKN = trimSelfConnecting(PKN, allTFs,MaintainInteractions)

PKN = MaxSubgraph(PKN)
allTFs = numpy.intersect1d(allTFs, PKN['target'])


PKN.to_csv('../preprocessing/preprocessed_data/PKN-Model_smaller.tsv', sep='\t', index=False)

