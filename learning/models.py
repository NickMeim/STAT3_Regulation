from __future__ import absolute_import, division
from tqdm.auto import tqdm
import torch
from torch.autograd import Variable, grad
import numpy as np
import pandas as pd
import sys
import random
import os

class Encoder(torch.nn.Module):
    def __init__(self, in_channel, hidden_layers, latent_dim,dropRate=0.1, activation=None, bias=True):

        super(Encoder, self).__init__()

        self.bias = bias
        self.num_hidden_layers = len(hidden_layers)
        self.bn = torch.nn.ModuleList()
        self.linear_layers = torch.nn.ModuleList()
        self.linear_layers.append(torch.nn.Linear(in_channel, hidden_layers[0], bias=bias))
        self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[0], momentum=0.6))
        for i in range(1, len(hidden_layers)):
            self.linear_layers.append(torch.nn.Linear(hidden_layers[i - 1], hidden_layers[i], bias=bias))
            self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[i], momentum=0.6))

        self.linear_latent_mu = torch.nn.Linear(hidden_layers[-1],
                                                latent_dim,
                                                bias=False)
        self.linear_latent_var = torch.nn.Linear(hidden_layers[-1],
                                                latent_dim,
                                                 bias=False)
        if activation is not None:
            self.activation = activation
        self.dropout = torch.nn.Dropout(dropRate)
        # self.drop_in = torch.nn.Dropout(0.5)

        self.N = torch.distributions.Normal(0, 1)
        self.N.loc = self.N.loc.cuda()  # hack to get sampling on the GPU
        self.N.scale = self.N.scale.cuda()

        self.init_emb()

    def init_emb(self):
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.fill_(0.0)


    def reparameterize(self, mu, log_var):

        if self.training:
            std = torch.exp(0.5 * log_var)
            eps = self.N.sample(std.shape)
            z = self.N.sample(mu.shape).mul(std).add_(mu)
            return z
        else:
            return mu

    def forward(self, x):
        # x = self.drop_in(x)
        for i in range(self.num_hidden_layers):
            x = self.linear_layers[i](x)
            x = self.bn[i](x)
            x = self.activation(x)
            x = self.dropout(x)
        z_mu = self.linear_latent_mu(x)
        z_log_var = self.linear_latent_var(x)
        z_latent = self.reparameterize(z_mu, z_log_var)

        return z_latent

    def L2Regularization(self, L2):

        weightLoss = 0.
        biasLoss = 0.
        for i in range(self.num_hidden_layers):
            weightLoss = weightLoss + L2 * torch.sum((self.linear_layers[i].weight)**2)
            if self.bias==True:
                biasLoss = biasLoss + L2 * torch.sum((self.linear_layers[i].bias)**2)
        L2Loss = biasLoss + weightLoss
        return(L2Loss)



class Decoder(torch.nn.Module):
    def __init__(self, latent_dim, hidden_layers, out_dim,dropRate=0.1, activation=None, bias=True):

        super(Decoder, self).__init__()

        self.bias = bias
        self.num_hidden_layers = len(hidden_layers)
        self.bn = torch.nn.ModuleList()
        self.linear_layers = torch.nn.ModuleList()
        self.linear_layers.append(torch.nn.Linear(latent_dim, hidden_layers[0], bias=bias))
        self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[0], momentum=0.6))
        for i in range(1, len(hidden_layers)):
            self.linear_layers.append(torch.nn.Linear(hidden_layers[i - 1], hidden_layers[i], bias=bias))
            self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[i], momentum=0.6))

        self.output_linear = torch.nn.Linear(hidden_layers[-1],
                                             out_dim,
                                             bias=False)

        if activation is not None:
            self.activation = activation
        self.dropout = torch.nn.Dropout(dropRate)

        self.init_emb()

    def init_emb(self):
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.fill_(0.0)

    def forward(self, x):
        for i in range(self.num_hidden_layers):
            x = self.linear_layers[i](x)
            x = self.bn[i](x)
            x = self.activation(x)
            x = self.dropout(x)

        output = self.output_linear(x)
        return output

    def L2Regularization(self, L2):

        weightLoss = 0.
        biasLoss = 0.
        for i in range(self.num_hidden_layers):
            weightLoss = weightLoss + L2 * torch.sum((self.linear_layers[i].weight)**2)
            if self.bias==True:
                biasLoss = biasLoss + L2 * torch.sum((self.linear_layers[i].bias)**2)
        L2Loss = biasLoss + weightLoss
        return(L2Loss)

class VAE(torch.nn.Module):

    def __init__(self, enc, dec):

        super(VAE, self).__init__()

        self.encoder = enc
        self.decoder = dec
        self.log_scale = torch.nn.Parameter(torch.Tensor([0.0]))

    def forward(self, x):
        z_latent = self.encoder(x)
        predicted = self.decoder(z_latent)

        return z_latent, predicted

    def encode(self, x):
        z_latent = self.encoder(x)
        return z_latent

    def decode(self, x):
        decoded_output = self.decoder(x)
        return decoded_output

    def L2Regularization(self, L2):

        encoderL2 = self.encoder.L2Regularization(L2)
        decoderL2 = self.decoder.L2Regularization(L2)

        L2Loss = encoderL2 + decoderL2
        return(L2Loss)