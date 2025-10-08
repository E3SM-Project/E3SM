import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class ANN_CNN(nn.Module):
    def __init__(self, idim, odim, hdim, stencil, dropout=0.0):
        super().__init__()

        self.idim = idim
        self.odim = odim
        self.hdim = hdim
        self.dropout_prob = dropout
        self.stencil = stencil
        self.fac = int(np.floor(0.5 * self.stencil))

        # assume normalized data as input
        # same activation for all layers
        # same dropout probabilities for all layers

        # Applying multiple 3x3 conv layers than just one stencilxstencil layer performs better
        if self.fac == 1:
            self.conv1 = nn.Conv2d(
                in_channels=idim, out_channels=idim, kernel_size=3, stride=1, padding=0
            )
            self.act_cnn = nn.ReLU()
            self.dropout0 = nn.Dropout(p=0.5 * self.dropout_prob)

        elif self.fac == 2:
            self.conv1 = nn.Conv2d(
                in_channels=idim, out_channels=idim, kernel_size=3, stride=1, padding=0
            )
            self.act_cnn = nn.ReLU()
            self.dropout0 = nn.Dropout(p=0.5 * self.dropout_prob)

            self.conv2 = nn.Conv2d(
                in_channels=idim, out_channels=idim, kernel_size=3, stride=1, padding=0
            )
            self.act_cnn2 = nn.ReLU()
            self.dropout0_2 = nn.Dropout(p=0.5 * self.dropout_prob)

        elif self.fac == 3:  # this is not ready yet, and might not be needed for my study
            self.conv1 = nn.Conv2d(
                in_channels=idim, out_channels=idim, kernel_size=self.stencil, stride=1, padding=0
            )
            self.act_cnn = nn.ReLU()
            self.dropout0 = nn.Dropout(p=0.5 * self.dropout_prob)

        # can define a block and divide it into blocks as well
        self.layer1 = nn.Linear(idim, hdim)  # ,dtype=torch.float16)
        self.act1 = nn.LeakyReLU()

        self.dropout = nn.Dropout(p=self.dropout_prob)

        self.layer2 = nn.Linear(hdim, hdim)
        self.act2 = nn.LeakyReLU()
        # -------------------------------------------------------
        self.layer3 = nn.Linear(hdim, hdim)
        self.act3 = nn.LeakyReLU()
        # -------------------------------------------------------
        self.layer4 = nn.Linear(hdim, hdim)
        self.act4 = nn.LeakyReLU()
        # --------------------------------------------------------
        self.layer5 = nn.Linear(hdim, hdim)
        self.act5 = nn.LeakyReLU()
        # -------------------------------------------------------
        self.layer6 = nn.Linear(hdim, 2 * odim)
        self.act6 = nn.LeakyReLU()

        self.output = nn.Linear(2 * odim, odim)

    def forward(self, x):
        if not torch.jit.is_scripting():
            # ignore this branch of code if torch is scripting
            if self.fac == 1:
                x = torch.squeeze(self.dropout0(self.act_cnn(self.conv1(x))))
            elif self.fac == 2:
                x = torch.squeeze(self.dropout0(self.act_cnn(self.conv1(x))))
                x = torch.squeeze(self.dropout0_2(self.act_cnn2(self.conv2(x))))

        # print(f'new shape: {x.shape}')
        x = self.dropout(self.act1(self.layer1(x)))
        x = self.dropout(self.act2(self.layer2(x)))
        x = self.dropout(self.act3(self.layer3(x)))
        x = self.dropout(self.act4(self.layer4(x)))
        x = self.dropout(self.act5(self.layer5(x)))
        x = self.dropout(self.act6(self.layer6(x)))
        x = self.output(x)

        return x

    def gaussian_dropout(self, x):
        S = x.shape
        vec = torch.normal(mean=torch.ones(S), std=torch.ones(S))
        return x * vec

    # calculates total number of learnable parameters
    def totalparams(self):
        param_size = 0
        for param in self.parameters():
            param_size += param.nelement()

        return param_size

    # computes total model size in MBs
    def totalsize(self):
        param_size = 0
        for param in self.parameters():
            param_size += param.nelement() * param.element_size()
        buffer_size = 0
        for buffer in self.buffers():
            buffer_size += buffer.nelement() * buffer.element_size()

        size_all_mb = (param_size + buffer_size) / 1024**2
        # print('model size: {:.3f}MB'.format(size_all_mb))

        return size_all_mb
