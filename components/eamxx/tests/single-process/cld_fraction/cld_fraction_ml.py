import os

import torch
import torch.nn as nn

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
    y21 = self.tot1(torch.cat((liq, y13_categorical), dim=1))
    y22 = self.relu(y21)
    y23 = self.tot2(y22)
    return y13_categorical, y23

model = None

def init ():
    global model

    # For this test, hard code nlevs, as well as pth file name/path
    nlevs = 72
    current_file_directory = os.path.dirname(os.path.abspath(__file__))
    model_file = f"{current_file_directory}/cldfrac_net_weights.pth"

    model = CldFracNet(nlevs,nlevs)
    model.load_state_dict(torch.load(model_file,map_location=torch.device('cpu')))

def main (ice_threshold, ice_4out_threshold,
          qi, liq_cld_frac,
          ice_cld_frac, tot_cld_frac,
          ice_cld_frac_4out, tot_cld_frac_4out):
    global model

    # Convert numpy inputs to torch arrays
    # Note: our pth model expects float32 arrays, so make sure we get the right dtype
    liq_pt = torch.tensor(liq_cld_frac, dtype=torch.float32)
    qi_pt  = torch.tensor(qi, dtype=torch.float32)

    # Set model in evaluation mode
    model.eval()

    with torch.no_grad(): # Disable gradient for inference
        # Run the emulator
        ice_out, tot_out = model(qi_pt,liq_pt)

        # Update inout numpy arrays inplace
        ice_cld_frac[:] = ice_out.cpu().numpy()
        tot_cld_frac[:] = tot_out.cpu().numpy()
