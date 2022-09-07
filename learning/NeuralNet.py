import torch
import numpy
import matplotlib.pyplot as plt
import bionetwork
import plotting
import pandas
import time
import copy
import argparse
import activationFunctions
from scipy import stats
from sklearn.preprocessing import StandardScaler

#parser = argparse.ArgumentParser(prog='KDs simulation')
#parser.add_argument('--leaveIn', action='store', default=None)
#args = parser.parse_args()
#curentId = int(args.leaveIn)
curentId = 0


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


class cellLayer(torch.nn.Module):
    def __init__(self, numberOfGenes):
        super(cellLayer, self).__init__()

        weights = torch.ones(numberOfGenes, requires_grad=True, dtype=torch.double)
        weights.data = 1e-3 * weights.data
        bias = torch.zeros(numberOfGenes, requires_grad=True, dtype=torch.double)
        self.weights = torch.nn.Parameter(weights)
        self.bias = torch.nn.Parameter(bias)

    def forward(self, dataCell):
        cellIn = dataCell * self.weights + self.bias

        # leaky cutoff, corresponds to leaky relu but in the oposit direction
        cellInFilter = cellIn.detach() > 0
        cellIn[cellInFilter] = 0.01 * cellIn[cellInFilter]
        return cellIn

    def signRegularization(self, MoAFactor):
        weightFilter = self.weights.detach() < 0
        return MoAFactor * torch.sum(torch.abs(self.weights[weightFilter]))

    def L2Regularization(self, L2):
        L2weight = torch.sum(torch.square(self.weights))
        L2bias = torch.sum(torch.square(self.bias))
        return L2 * (L2weight + L2bias)

class fullModel(torch.nn.Module):
    def __init__(self, nodeNames, inNames, outName, networkList, modeOfAction, inputAmplitude,projectionAmplitude):
        super(fullModel, self).__init__()

        bionetParams = bionetwork.trainingParameters(iterations=150, clipping=1, targetPrecision=1e-6, leak=0.01)

        self.cellModel = cellLayer(len(nodeNames))
        self.inputLayer = bionetwork.projectInput(nodeNames, inNames, inputAmplitude, torch.double)
        self.signalingModel = bionetwork.bionet(networkList, len(nodeNames), modeOfAction, bionetParams, 'MML',
                                                torch.double)
        self.projectionLayer = bionetwork.projectOutput(nodeNames, outName, projectionAmplitude, torch.double)

    def forward(self, dataIn,dataCell, noiseLevel=0):
        Yin = self.inputLayer(dataIn) + self.cellModel(dataCell)
        Yin = Yin + noiseLevel * torch.randn(Yin.shape)
        YhatFull = self.signalingModel(Yin)
        Yhat = self.projectionLayer(YhatFull)

        return Yhat, YhatFull

    def L2Regularization(self, L2, lbRegularization=0.001):
        cellL2 = self.cellModel.L2Regularization(L2)

        # Signaling model L2
        absFilter = torch.abs(self.signalingModel.weights.detach()) > lbRegularization
        weightLoss = L2 * torch.sum(torch.square(self.signalingModel.weights[absFilter.detach()]))
        biasLoss = L2 * torch.sum(torch.square(self.signalingModel.bias))

        L2Loss = weightLoss + biasLoss + cellL2
        return L2Loss

    def signRegularization(self, MoAFactor):
        cellSign = self.cellModel.signRegularization(MoAFactor)
        signalingSign = self.signalingModel.signRegularization(MoAFactor)
        signConstraints = cellSign + signalingSign
        return signConstraints


inputAmplitude = 1.2
projectionAmplitude = 0.1

# Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/KEGG-Model.tsv')

# Load input output data
#Load samples
train_samples = pandas.read_csv('../data/10fold_cross_validation/KEGG/train_sample_' + str(curentId) + '.csv', low_memory=False,index_col=0)
train_samples = train_samples.sample(frac=1).reset_index(drop=True)
#Load TFs
TFsData = pandas.read_csv('../results/trimmed_shrnas_tf_activities_kegg.tsv',sep='\t', low_memory=False,index_col=0)
trainTFs = TFsData.loc[train_samples.sig_id.values,:]
outName = list(trainTFs.columns)
#Load KDs
allKds = pandas.read_csv('../preprocessing/preprocessed_data/all_filtered_Kds_kegg.tsv',sep='\t', low_memory=False,index_col=0)
allKds = allKds.set_index('sig_id')
inName = list(allKds.columns)
trainKds = allKds.loc[train_samples.sig_id.values,:]
# trainKds = trainKds.T # only for cell-line specific models
# missingValues = numpy.setdiff1d(nodeNames, trainKds.index.values) # only for cell-line specific models
# df = pandas.DataFrame(numpy.zeros((len(missingValues), trainKds.shape[1])), index=missingValues, columns=trainKds.columns)# only for cell-line specific models
# trainKds = trainKds.append(df)# only for cell-line specific models
# trainKds = trainKds.loc[nodeNames,:] # only for cell-line specific models
# trainKds = trainKds.T # only for cell-line specific models
# inName = list(trainKds.columns)  # only for cell-line specific models
#Load cell line data
cellLineMember = pandas.read_csv('../preprocessing/preprocessed_data/all_filtered_cells_kegg.tsv', sep='\t', low_memory=False, index_col=0)
cellLineMember = cellLineMember.set_index('sig_id')
TrainCellLineMember = cellLineMember.loc[train_samples.sig_id.values,:]
cellLineLevels = pandas.read_csv('../data/CCLE/trimmed_ccle_kegg.tsv', sep='\t', low_memory=False, index_col=0)
# cellLineLevels = cellLineLevels.T # only for cell-line specific models
# TrainCellLineMember = TrainCellLineMember.loc[:,'HT29'] # only for cell-line specific models
# missingValues = numpy.setdiff1d(nodeNames, cellLineLevels.index.values) # only for cell-line specific models
# df = pandas.DataFrame(numpy.zeros((len(missingValues), cellLineLevels.shape[1])), index=missingValues, columns=cellLineLevels.columns)# only for cell-line specific models
# cellLineLevels = cellLineLevels.append(df)# only for cell-line specific models
# cellLineLevels = cellLineLevels.loc[nodeNames,:]# only for cell-line specific models
# cellLineLevels = cellLineLevels.loc[:,'HT29'] # only for cell-line specific models
# scaler = StandardScaler()
# cellLineLevels_scaled = pandas.DataFrame(scaler.fit_transform(cellLineLevels))
# cellLineLevels_scaled.index = cellLineLevels.index
# cellLineLevels_scaled.columns = cellLineLevels.columns
cellLineLevels_scaled = cellLineLevels.T
missingValues = numpy.setdiff1d(nodeNames, cellLineLevels_scaled.index.values)
#Zero padding:
df = pandas.DataFrame(numpy.zeros((len(missingValues), cellLineLevels_scaled.shape[1])), index=missingValues, columns=cellLineLevels_scaled.columns)
cellLineLevels_scaled = cellLineLevels_scaled.append(df)
cellLineLevels_scaled = cellLineLevels_scaled.loc[nodeNames,:]
#cellLineLevels_scaled = cellLineLevels_scaled/cellLineLevels_scaled.max() # scale between zero and one
geneData = cellLineLevels_scaled.values.dot(TrainCellLineMember.values.T).T
# cellLineLevels = cellLineLevels/cellLineLevels.max()# only for cell-line specific models
# geneData = cellLineLevels.values.reshape(cellLineLevels.values.shape[0],1).dot(TrainCellLineMember.values.reshape(1,TrainCellLineMember.shape[0])).T # only for cell-line specific models


# Build model
model = fullModel(nodeNames, inName, outName, networkList, modeOfAction, inputAmplitude, projectionAmplitude)
model.signalingModel.preScaleWeights()
model.inputLayer.weights.requires_grad = False
#model.projectionLayer.weights.requires_grad = False
model.signalingModel.bias.data[numpy.isin(nodeNames, outName)] = 0.5 # put all TFs in 0.5 to begin

criterion = torch.nn.MSELoss()

sampleName = trainKds.index.values
X = torch.tensor(trainKds.values, dtype=torch.double)# *0.5 only for cell-line specific models
#X = torch.tensor(trainKds.values, dtype=torch.double) + torch.tensor(geneData, dtype=torch.double) # only for cell-line specific models
Xcell = torch.tensor(geneData, dtype=torch.double)
# Xcell = torch.tensor(geneData, dtype=torch.double).detach() # only for cell-line specific models
Y = torch.tensor(trainTFs.values.copy(), dtype=torch.double)

# Setup optimizer
MoAFactor = 0.1
spectralFactor = 1e-3
noiseLevel = 1e-3
batchSize = 230
maxIter = 20000
L2 = 1e-6

referenceState = copy.deepcopy(model.state_dict())
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=0)
#scheduler = torch.optim.lr_scheduler.StepLR(optimizer,
#                                            step_size=5000,
#                                            gamma=0.5)
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
    curLr = bionetwork.oneCycle(e, maxIter, maxHeight=2e-3, startHeight=5e-4, endHeight=5e-3, peak=2000)
    optimizer.param_groups[0]['lr'] = curLr
    curLoss = []
    curEig = []
    trainloader = bionetwork.getSamples(N, batchSize)
    for dataIndex in trainloader:
        model.train()
        model.signalingModel.weights.data = model.signalingModel.weights.data + 1e-8 * torch.randn(model.signalingModel.weights.shape)  # breaks potential symmetries
        optimizer.zero_grad()
        
        dataIn = X[dataIndex, :].view(len(dataIndex), X.shape[1])
        dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])
        dataCell = Xcell[dataIndex, :].view(len(dataIndex), Xcell.shape[1])

        #dataCell = dataCell + noiseLevel * torch.randn(dataCell.shape)
        Yhat, YhatFull = model(dataIn, dataCell, noiseLevel)
        # Yhat, YhatFull = model(dataIn, noiseLevel)

        curState[dataIndex, :] = YhatFull.detach()

        fitLoss = criterion(dataOut, Yhat)

        stateLoss = 1e-5 * uniformLoss(curState, dataIndex, YhatFull, maxConstraintFactor=10.)
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
    #scheduler.step()

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
    if e % 1000 == 0:
        print(outString)

    # if numpy.logical_and(e % 250 == 0, e>0):
    #    optimizer.state = resetState.copy()

model.eval()
Yhat, YhatFull = model(X,Xcell)
# Yhat, YhatFull = model(X)
torch.save(model, '../results/crossValRes/eptest/model_' + str(curentId) + '.pt')

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
plt.savefig('../figures/crossValRes/eptest/training_loss_' + str(curentId) + '.png', dpi=600)

plt.figure()
plt.plot(T, mseLoss)
plt.plot([0, len(T)], numpy.array([1, 1]) * mLoss.item(), 'black', linestyle='--')
plt.xlim([0, len(T)])
plt.ylabel('MSE')
plt.xlabel('Epoch')
plt.yscale('log')
plt.savefig('../figures/crossValRes/eptest/training_mse_' + str(curentId) + '.png', dpi=600)

plt.figure()
plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
plt.savefig('../figures/crossValRes/eptest/fit_performance_' + str(curentId) + '.png', dpi=600)
plt.rcParams["figure.figsize"] = (12, 10)
plt.figure()
rank = plotting.compareAllTFs(Yhat, Y, outName)
plt.savefig('../figures/crossValRes/eptest/perTF_fit_train_' + str(curentId) + '.png', dpi=600)
# plt.figure()
# rank = plotting.compareAllTFs(Yhat.T, Y.T, sampleName)
#
# plt.savefig('../figures/crossValRes/eptest/perSample_fit_train_' + str(curentId) + '.png', dpi=600)
# plotting.displayData(Y, sampleName, outName)
# plt.savefig('../figures/crossValRes/heatmap_' + str(curentId) + '.png', dpi=600)

fitData = pandas.DataFrame(Yhat.detach().numpy(), index=sampleName, columns=outName)
fitData.to_csv('../results/crossValRes/eptest/fit_' + str(curentId) + '.tsv', sep='\t')

plt.figure()
plt.hist(spectral_r)
plt.savefig('../figures/crossValRes/eptest/spectral_r_' + str(curentId) + '.png', dpi=600)
## Validation
val_samples = pandas.read_csv('../data/10fold_cross_validation/KEGG/val_sample_' + str(curentId) + '.csv', low_memory=False,index_col=0)
#Load TFs
valTFs = TFsData.loc[val_samples.sig_id.values,:]
#Load KDs
valKds = allKds.loc[val_samples.sig_id.values,:]
# valKds = valKds.T # only for cell-line specific models
# missingValues = numpy.setdiff1d(nodeNames, valKds.index.values) # only for cell-line specific models
# df = pandas.DataFrame(numpy.zeros((len(missingValues), valKds.shape[1])), index=missingValues, columns=valKds.columns)# only for cell-line specific models
# valKds = valKds.append(df)# only for cell-line specific models
# valKds = valKds.loc[nodeNames,:] # only for cell-line specific models
# valKds = valKds.T # only for cell-line specific models
# inName = list(valKds.columns)  # only for cell-line specific models
#Load cell line data
ValCellLineMember = cellLineMember.loc[val_samples.sig_id.values,:]
geneData = cellLineLevels_scaled.values.dot(ValCellLineMember.values.T).T

sampleName = val_samples.sig_id.values
X = torch.tensor(valKds.values, dtype=torch.double)
Xcell = torch.tensor(geneData, dtype=torch.double)
Y = torch.tensor(valTFs.values.copy(), dtype=torch.double)

Yhat, YhatFull = model(X,Xcell)
# Yhat, YhatFull = model(X)
plt.figure()
plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
plt.savefig('../figures/crossValRes/eptest/validation_performance_' + str(curentId) + '.png', dpi=600)

plt.figure()
rank = plotting.compareAllTFs(Yhat, Y, outName)
plt.savefig('../figures/crossValRes/eptest/perTF_fit_validation_' + str(curentId) + '.png', dpi=600)
# plt.figure()
# rank = plotting.compareAllTFs(Yhat.T, Y.T, sampleName)
# plt.savefig('../figures/crossValRes/eptest/perSample_fit_validation_' + str(curentId) + '.png', dpi=600)

