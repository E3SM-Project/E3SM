#!/usr/bin/env python3
import os

import torch
import torch.nn as nn
import torch_mlir
from torch_mlir.compiler_utils import TensorPlaceholder
from torch_mlir import torchscript
from lapis import KokkosBackend

class CldFracNet(nn.Module):
  def __init__(self, input_size, output_size, neuron_count=64):
    super(CldFracNet, self).__init__()
    # emulate cld_ice = (qi > 1e-5)
    self.ice1 = nn.Linear(input_size, neuron_count)
    self.ice2 = nn.Linear(neuron_count, output_size)
    # emulate cld_tot = max(cld_ice, cld_liq)
    self.tot1 = nn.Linear(input_size*2, neuron_count)
    self.tot2 = nn.Linear(neuron_count, output_size)
    # a relu for fun
    self.relu = nn.ReLU()
    # sigmoid for categorical ice output
    self.sigmoid = nn.Sigmoid()

  def forward(self, qi, liq):
    # First, compute cld_ice from qi
    y11 = self.ice1(qi)
    y12 = self.relu(y11)
    y13 = self.ice2(y12)
    # Apply sigmoid to get probabilities
    y13_probabilities = self.sigmoid(y13)

    # During training, use straight-through estimator for gradients
    # During inference, use hard binary values
    if self.training:
        # Straight-through estimator: forward pass uses binary, backward pass uses sigmoid
        y13_binary = (y13_probabilities > 0.5).float()
        y13_categorical = y13_binary - y13_probabilities.detach() + y13_probabilities
    else:
        # During inference, use hard binary values
        y13_categorical = (y13_probabilities > 0.5).float()

    # Now compute cld_tot from cld_ice and cld_liq
    y21 = self.tot1(torch.cat((liq, y13_categorical), dim=0))
    y22 = self.relu(y21)
    y23 = self.tot2(y22)
    return y13_categorical, y23

def main ():

    # For this test, hard code nlevs, as well as pth file name/path
    nlevs = 72
    current_file_directory = os.path.dirname(os.path.abspath(__file__))
    model_file = f"{current_file_directory}/cldfrac_net_weights.pth"

    model = CldFracNet(nlevs,nlevs)
    model.load_state_dict(torch.load(model_file,map_location=torch.device('cpu')))

    qi = torch.ones((nlevs))
    liq = torch.ones((nlevs))
    mlir_module = torchscript.compile(model, (qi, liq), output_type='linalg-on-tensors')
    with open("cldfrac.mlir",'w') as fd:
        fd.write(str(mlir_module))

if __name__ == "__main__":
    main()

