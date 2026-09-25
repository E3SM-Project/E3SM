"""A torch.nn model the Python backend test runs. Needs torch."""

import torch


def create_emulator(config):
    return Emulator(config)


class Net(torch.nn.Module):
    """ReLU(W x + b) over the last dimension, 3 features in and 2 out.

    The weights are fixed and small, so results are exact in float32.
    """

    def __init__(self):
        super().__init__()
        self.linear = torch.nn.Linear(3, 2)
        with torch.no_grad():
            self.linear.weight.copy_(
                torch.tensor([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            )
            self.linear.bias.copy_(torch.tensor([0.5, -0.5]))

    def forward(self, x):
        return torch.relu(self.linear(x))


class Emulator:
    def __init__(self, config):
        device = config.get("device", "auto")
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device)
        self.net = Net().to(self.device).eval()
        if config["verbose"]:
            print(f"[python_torch_fixture] running on {self.device}", flush=True)

    def infer(self, inputs, outputs):
        # float64 model memory -> the network's float32, on its device
        x = torch.tensor(inputs["x"], dtype=torch.float32, device=self.device)
        with torch.no_grad():
            y = self.net(x)
        outputs["y"][:] = y.double().cpu().numpy()
