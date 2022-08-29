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

# Probably large enough batches to use batch level regularization
# def uniformLoss(data, targetMin = 0, targetMax = 0.99, maxConstraintFactor = 1):
#     targetMean = (targetMax-targetMin)/2
#     targetVar= (targetMax-targetMin)**2/12
#
#     nodeMean = torch.mean(data, dim=0)
#     nodeVar = torch.mean(torch.square(data-nodeMean), dim=0)
#     maxVal, _ = torch.max(data, dim=0)
#     minVal, _ = torch.min(data, dim=0)
#
#     meanLoss = torch.sum(torch.square(nodeMean - targetMean))
#     varLoss =  torch.sum(torch.square(nodeVar - targetVar))
#     maxLoss = torch.sum(torch.square(maxVal - targetMax))
#     minloss = torch.sum(torch.square(minVal- targetMin))
#     maxConstraint = -maxConstraintFactor * torch.sum(maxVal[maxVal.detach()<=0]) #max value should never be negative
#
#     loss = meanLoss + varLoss + minloss + maxLoss + maxConstraint
#     return loss

def uniformLoss(curState, dataIndex, YhatFull, targetMin = 0, targetMax = 0.99, maxConstraintFactor = 10):
    data = curState.detach().clone()
    data[dataIndex, :] = YhatFull

    targetMean = (targetMax-targetMin)/2
    targetVar= (targetMax-targetMin)**2/12

    factor = 1
    meanFactor = factor
    varFactor = factor
    minFactor = factor
    maxFactor = factor
    maxConstraintFactor = factor * maxConstraintFactor

    nodeMean = torch.mean(data, dim=0)
    nodeVar = torch.mean(torch.square(data-nodeMean), dim=0)
    maxVal, _ = torch.max(data, dim=0)
    minVal, _ = torch.min(data, dim=0)

    meanLoss = meanFactor * torch.sum(torch.square(nodeMean - targetMean))
    varLoss =  varFactor * torch.sum(torch.square(nodeVar - targetVar))
    maxLoss = maxFactor * torch.sum(torch.square(maxVal - targetMax))
    minloss = minFactor * torch.sum(torch.square(minVal- targetMin))
    maxConstraint = -maxConstraintFactor * torch.sum(maxVal[maxVal.detach()<=0]) #max value should never be negative

    loss = meanLoss + varLoss + minloss + maxLoss + maxConstraint
    return loss

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

bionetParams = bionetwork.trainingParameters(iterations = 150, clipping=1, targetPrecision=1e-6, leak=0.01)

class fullModel(torch.nn.Module):
    def __init__(self, nodeNames, outName, networkList, modeOfAction,projectionAmplitude):
        super(fullModel, self).__init__()

        bionetParams = bionetwork.trainingParameters(iterations=150, clipping=1, targetPrecision=1e-6, leak=0.01)

        self.signalingModel = bionetwork.bionet(networkList, len(nodeNames), modeOfAction, bionetParams, 'MML',
                                                torch.double)
        self.projectionLayer = bionetwork.projectOutput(nodeNames, outName, projectionAmplitude, torch.double)

    def forward(self, dataIn, noiseLevel=0):
        Yin = dataIn + noiseLevel * torch.randn(dataIn.shape)
        YhatFull = self.signalingModel(Yin)
        Yhat = self.projectionLayer(YhatFull)

        return Yhat, YhatFull

    def L2Regularization(self, L2, lbRegularization=0.001):

        # Signaling model L2
        absFilter = torch.abs(self.signalingModel.weights.detach()) > lbRegularization
        weightLoss = L2 * torch.sum(torch.square(self.signalingModel.weights[absFilter.detach()]))
        biasLoss = L2 * torch.sum(torch.square(self.signalingModel.bias))

        L2Loss = weightLoss + biasLoss
        return L2Loss

    def signRegularization(self, MoAFactor):
        signalingSign = self.signalingModel.signRegularization(MoAFactor)
        return signalingSign


#inputAmplitude = 1.2
projectionAmplitude = 0.1

# Load data
X = torch.load('inputX_0.pt')
Y = torch.load('KORef_0.pt')
Yfull = torch.load('KOFullRef_0.pt')

Xtest = torch.load('inputX_1.pt')
Ytest = torch.load('KORef_1.pt')
Yfulltest = torch.load('KOFullRef_1.pt')

# Build model
model = fullModel(nodeNames, outName, networkList, modeOfAction, projectionAmplitude)
model.signalingModel.preScaleWeights()
model.signalingModel.bias.data[numpy.isin(nodeNames, outName)] = 0.5 # put all TFs in 0.5 to begin
# model.signalingModel.bias.data[numpy.isin(nodeNames, inName)] = 1  #Begin with signalal from all ligands

criterion = torch.nn.MSELoss()

# Setup optimizer
MoAFactor = 0.1
spectralFactor = 1e-3
noiseLevel = 1e-3
batchSize = 25
maxIter = 1000
L2 = 1e-5

referenceState = copy.deepcopy(model.state_dict())
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=0)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer,
                                           step_size=200,
                                           gamma=0.5)
resetState = optimizer.state.copy()

mLoss = criterion(torch.mean(Y, dim=0) * torch.ones(Y.shape), Y)
print(mLoss)

N = X.shape[0]

e = 0
trainLoss = []
mseLoss = []
spectral_r = []
curState = torch.rand((Y.shape[0], model.signalingModel.bias.shape[0]), dtype=torch.double, requires_grad=False)

for e in range(e, maxIter):
    curLoss = []
    curEig = []
    trainloader = bionetwork.getSamples(N, batchSize)
    for dataIndex in trainloader:
        model.train()
        model.signalingModel.weights.data = model.signalingModel.weights.data + 1e-8 * torch.randn(model.signalingModel.weights.shape)  # breaks potential symmetries
        optimizer.zero_grad()
        
        dataIn = X[dataIndex, :].view(len(dataIndex), X.shape[1])
        dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])

        Yhat, YhatFull = model(dataIn, noiseLevel)

        curState[dataIndex, :] = YhatFull.detach()

        fitLoss = criterion(dataOut, Yhat)

        stateLoss = 1e-5 * uniformLoss(curState, dataIndex, YhatFull,maxConstraintFactor=50.)
        #stateLoss = 1e-6 * uniformLoss(YhatFull)
        L2Loss = model.L2Regularization(L2)
        signConstraint = model.signRegularization(MoAFactor)

        spectralRadiusLoss, spectralRadius = bionetwork.spectralLoss(model.signalingModel, YhatFull.detach(),
                                                                     model.signalingModel.weights, expFactor=5)

        projectionLoss = 1e-6 * torch.sum(torch.square(model.projectionLayer.weights - projectionAmplitude))

        loss = fitLoss + signConstraint + stateLoss + L2Loss  + spectralFactor * spectralRadiusLoss + projectionLoss

        loss.backward()

        optimizer.step()

        curEig.append(spectralRadius.item())
    spectral_r.append(numpy.max(curEig))
    scheduler.step()

    # stats = plotting.storeProgress(stats, e, curLoss, curEig, curLr, violations=model.signalingModel.getNumberOfViolations())

    outString = 'Epoch={:.0f}/{:.0f}'.format(e + 1, maxIter)
    outString += ',Loss={:.4f}'.format(loss.item())
    outString += ',MSE={:.4f}'.format(fitLoss.item())
    outString += ',State Loss={:.4f}'.format(stateLoss.item())
    outString += ',L2 Loss={:.5f}'.format(L2Loss.item())
    outString += ',Project L2={:.4f}'.format(projectionLoss.item())
    outString += ',Sign Loss={:.4f}'.format(signConstraint.item())
    outString += ',Spectral Loss={:.4f}'.format(spectralRadiusLoss.item())
    trainLoss.append(loss.item())
    mseLoss.append(fitLoss.item())
    if e % 50 == 0:
        print(outString)

    if numpy.logical_and(e % 250 == 0, e>0):
       optimizer.state = resetState.copy()

model.eval()
Yhat, YhatFull = model(X)
torch.save(model, 'synthModel_control.pt')

# %%
spectralCapacity = numpy.exp(numpy.log(1e-6) / model.signalingModel.param['iterations'])

plt.rcParams["figure.figsize"] = (9, 6)

T = numpy.array(range(len(trainLoss)))

plt.figure()
plt.plot(T, trainLoss)
plt.xlim([0, len(T)])
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.yscale('log')
#plt.show()
plt.savefig('training_loss_0.png', dpi=600)

plt.figure()
plt.plot(T, mseLoss)
plt.plot([0, len(T)], numpy.array([1, 1]) * mLoss.item(), 'black', linestyle='--')
plt.xlim([0, len(T)])
plt.ylabel('MSE')
plt.xlabel('Epoch')
plt.yscale('log')
#plt.show()
plt.savefig('training_mse_0.png', dpi=600)

plt.figure()
plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
#plt.show()
plt.savefig('fit_performance_0.png', dpi=600)

plt.rcParams["figure.figsize"] = (12, 10)
plt.figure()
rank = plotting.compareAllTFs(Yhat, Y, outName)
#plt.show()
plt.savefig('perTF_fit_train_0.png', dpi=600)
plt.figure()
rank = plotting.compareAllTFs(Yhat.T, Y.T, sampleName)
#plt.show()
plt.savefig('perSample_fit_train_0.png', dpi=600)
# plotting.displayData(Y, sampleName, outName)
# plt.show()
# plt.savefig('heatmap.png', dpi=600)

fitData = pandas.DataFrame(Yhat.detach().numpy(), index=sampleName, columns=outName)
fitData.to_csv('fit_0.tsv', sep='\t')

plt.figure()
plt.hist(spectral_r)
plt.savefig('spectral_r_0.png', dpi=600)

# ## Validation
# #### Needs to implemet that
Yhat, YhatFull = model(Xtest)
plt.figure()
plt.scatter(Yhat.detach().numpy(), Ytest.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Ytest)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
# plt.show()
plt.savefig('validation_performance.png', dpi=600)
#
plt.figure()
rank = plotting.compareAllTFs(Yhat, Ytest, outName)
plt.savefig('perTF_fit_validation.png', dpi=600)
plt.figure()
rank = plotting.compareAllTFs(Yhat.T, Ytest.T, sampleName)
plt.savefig('perSample_fit_validation.png', dpi=600)

