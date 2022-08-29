import torch
import numpy
import matplotlib.pyplot as plt
import bionetwork
import pandas
from scipy.stats import pearsonr
import seaborn as sns
import copy
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib.colors import ListedColormap
import plotting
import argparse
#parser = argparse.ArgumentParser(prog='KDs simulation')
#parser.add_argument('--leaveIn', action='store', default=None)
#args = parser.parse_args()
#curentId = int(args.leaveIn)
curentId = 0


def loadModel(refModel, fileName):
    # work around to copy weights and bias values
    # because class structure has been updated since the run
    curModel = torch.load(fileName)
    model = copy.deepcopy(refModel)
    model.network.weights.data = curModel.network.weights.data.clone()
    model.network.bias.data = curModel.network.bias.data.clone()
    return model


def hashArray(array, hashMap):
    outputArray = array.copy()
    for i in range(len(array)):
        outputArray[i] = hashMap[array[i]]
    return outputArray


N = 2
simultaniousInput = 5
inputAmplitude = 3
projectionAmplitude = 1.2

# Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/KEGGnet-Model.tsv')
annotation = pandas.read_csv('data/KEGGnet-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))

# Subset input and output to intersecting nodes
inName = annotation.loc[annotation['ligand'], 'code'].values
outName = annotation.loc[annotation['TF'], 'code'].values
inName = numpy.intersect1d(nodeNames, inName)
outName = numpy.intersect1d(nodeNames, outName)
inNameGene = [uniprot2gene[x] for x in inName]
outNameGene = [uniprot2gene[x] for x in outName]
internalNodes = numpy.logical_not(numpy.logical_or(numpy.isin(nodeNames, inName), numpy.isin(nodeNames, outName)))
nodeNameGene = [uniprot2gene[x] for x in nodeNames]
sampleName = ['no_KO'] + list(numpy.array(nodeNameGene)[numpy.where(internalNodes)])

bionetParams = bionetwork.trainingParameters(iterations = 150, clipping=1, leak=0.01)
parameterizedModel = bionetwork.model(networkList, nodeNames, modeOfAction, inputAmplitude, projectionAmplitude, inName, outName, bionetParams)
parameterizedModel = bionetwork.loadParam('equationParams.txt', parameterizedModel, nodeNames)
parameterizedModel.projectionLayer.weights.data = parameterizedModel.projectionLayer.weights.data/projectionAmplitude

# Generate test data
Xtest = torch.zeros(N, len(inName), dtype=torch.double)
for i in range(1,N):  # skip 0 to include a ctrl sample i.e. zero input
    Xtest[i, numpy.random.randint(0, len(inName), simultaniousInput)] = torch.rand(simultaniousInput,
                                                                                   dtype=torch.double)
Ytest, YfullRef = parameterizedModel(Xtest)
Ytest = Ytest.detach()

#%%
#Knock out test
plt.rcParams["figure.figsize"] = (5, 5)
knockOutLevel = -5

nrOfKO = N
nrInternalNodes = sum(internalNodes)

for i in range(nrOfKO):
    inputX = parameterizedModel.inputLayer(Xtest[i,:].reshape(1,-1)).detach()
    inputX = inputX.repeat(sum(internalNodes)+1, 1)
    inputX[1:inputX.shape[0], :] = inputX[1:inputX.shape[0], :]  + knockOutLevel * torch.eye(len(nodeNames))[internalNodes,:]
    KOFullRef = parameterizedModel.network(inputX)
    KORef = parameterizedModel.projectionLayer(KOFullRef).detach()
    KOFullRef = KOFullRef.detach()
    KORef[KORef<0] = 0.
    KORef[KORef > 1] = 1.
    KOFullRef[KOFullRef < 0] = 0.
    KOFullRef[KOFullRef > 1] = 1.
    torch.save(KOFullRef,'KOFullRef_'+str(i)+'.pt')
    torch.save(KORef, 'KORef_' + str(i) + '.pt')
    torch.save(inputX.detach(), 'inputX_' + str(i) + '.pt')
    plotting.displayData(KORef, sampleName, outNameGene)
    plt.savefig('KOsignal_'+str(i)+'.png', dpi=600)
    print(KORef.max().item())
    print(KORef.min().item())
plt.figure()
plt.hist(KORef.flatten().numpy())
plt.show()
plt.figure()
plt.hist(KOFullRef.flatten().numpy())
plt.show()
