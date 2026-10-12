"""SDM-trained warm-rain emulator (ERF super-droplet LES), for P3's warm-rain emulator.

The model file (sdm_warm_rain.pt) holds the weights and everything specific to the
model: normalization, gates, training envelope and number-rate constants. A
retrained model is a file swap.

Used with warm_rain_emulator_python_module: sdm_warm_rain, and
warm_rain_emulator_model_path: the model file.

Contract (see eamxx_p3_warm_rain_emulator.cpp); arrays are (ncol, nlev) doubles:
  inputs : qc, nc, qr, nr [kg/kg, #/kg; cell averages, dry], rho [kg/m3, dry]
  outputs: qc2qr_autoconv_tend, qc2qr_accret_tend, ncautr, nc2nr_autoconv_tend,
           nc_accret_tend, nc_selfcollect_tend, nr_selfcollect_tend (cell averages,
           with P3's signs), and use_cloud, use_rain (1 inside the training envelope)
"""
import math
import warnings
from typing import Dict, Tuple

import numpy as np
import torch
from torch import nn

SUPPORTED_FORMAT = 1
INPUTS = ('qc', 'nc', 'qr', 'nr', 'rho')
OUTPUTS = ('qc2qr_autoconv_tend', 'qc2qr_accret_tend', 'ncautr', 'nc2nr_autoconv_tend',
           'nc_accret_tend', 'nc_selfcollect_tend', 'nr_selfcollect_tend', 'use_cloud', 'use_rain')
_ACT = dict(tanh=nn.Tanh, relu=nn.ReLU, silu=nn.SiLU, softplus=nn.Softplus)


class SdmWarmRain(nn.Module):
    """compact: run the network only on the cells where P3 uses its rates (inside the
    training envelope), and return zeros elsewhere. Most cells are clear or outside
    the envelope, and P3 ignores the rates there, so this only saves time. With
    compact=False every cell is computed, as in the model's reference vectors."""

    def __init__(self, model_file: str, compact: bool = True):
        super().__init__()
        self.compact = compact
        blob = torch.load(model_file, map_location='cpu', weights_only=True)
        cfg = blob['config']
        if cfg['format_version'] != SUPPORTED_FORMAT:
            raise RuntimeError(f"sdm_warm_rain: model format {cfg['format_version']} is not supported")
        names = [n for n, _ in cfg['contract']['inputs']]
        if names != ['qc', 'Nc', 'qr', 'Nr'] or \
           cfg['transforms'] != dict(input='log10_floor_standardize', output='expm1_logstd_scale'):
            raise RuntimeError(f"sdm_warm_rain: model contract {names} / {cfg['transforms']} is not supported")

        a = cfg['architecture']
        layers, n = [], a['n_in']
        for w in a['hidden']:
            layers += [nn.Linear(n, w), _ACT[a['activation']]()]
            n = w
        layers += [nn.Linear(n, a['n_out']), _ACT[a['output_activation']]()]
        self.net = nn.Sequential(*layers)
        for name in ('x_mean', 'x_std', 'floors', 'y_scale', 'y_log_std'):
            self.register_buffer(name, torch.zeros_like(blob['state_dict'][name]))
        self.load_state_dict(blob['state_dict'])
        self.requires_grad_(False)
        self.eval()

        g, c, r, nr = cfg['gates'], cfg['envelope']['cloud'], cfg['envelope']['rain'], cfg['number_rates']
        self.qc_gt, self.qr_gt = float(g['qc_gt']), float(g['qr_gt'])
        self.c_qc_max, self.c_nc_min, self.c_nc_max = float(c['qc_max']), float(c['nc_min']), float(c['nc_max'])
        self.c_qr_min, self.c_qr_max, self.c_nr_max = float(c['qr_min']), float(c['qr_max']), float(c['nr_max'])
        self.r_qr_max, self.r_nr_max = float(r['qr_max']), float(r['nr_max'])
        self.m_star = 4.0 / 3.0 * math.pi * 1000.0 * float(nr['embryo_radius_m']) ** 3
        self.drops_per_embryo = float(nr['cloud_drops_per_embryo'])
        self.ac_n_factor = float(nr['ac_n_factor'])
        # The rates are means over this period [s]: the model is meant for time steps close to it
        self.target_period_s = float(cfg['contract']['target_period_s'])

    def forward(self, qc: torch.Tensor, nc: torch.Tensor, qr: torch.Tensor, nr: torch.Tensor,
                rho: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor,
                                            torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor,
                                            torch.Tensor]:
        shape = qc.shape
        qc, nc, qr, nr, rho = [v.reshape(-1).double() for v in (qc, nc, qr, nr, rho)]
        qc_v, nc_v, qr_v, nr_v = qc * rho, nc * rho, qr * rho, nr * rho

        use_cloud = (qc_v > self.qc_gt) & (qc_v <= self.c_qc_max) & (nc_v >= self.c_nc_min) & \
                    (nc_v <= self.c_nc_max) & (qr_v > self.c_qr_min) & (qr_v <= self.c_qr_max) & \
                    (nr_v <= self.c_nr_max)
        use_rain = (qr_v <= self.r_qr_max) & (nr_v <= self.r_nr_max)

        # The cells to compute: all of them, or those where P3 uses a rate
        # (rain self-collection is zero below its gate anyway)
        if self.compact:
            cells = torch.nonzero(use_cloud | (use_rain & (qr_v > self.qr_gt))).squeeze(1)
        else:
            cells = torch.arange(qc.numel(), device=qc.device)

        # The network works in single precision, on per-volume quantities
        x = torch.stack([qc_v[cells], nc_v[cells], qr_v[cells], nr_v[cells]], 1).float()
        z = (torch.log10(torch.maximum(x, self.floors)) - self.x_mean) / self.x_std
        p = torch.expm1(self.net(z) * self.y_log_std) * self.y_scale
        cloud = x[:, 0] > self.qc_gt
        gates = torch.stack([cloud, cloud & (x[:, 2] > 0.0), cloud, x[:, 2] > self.qr_gt], 1)
        raw = (p * gates.to(p.dtype)).double() / rho[cells][:, None]

        qcc, ncc = qc[cells], nc[cells]
        nc_over_qc = torch.where(qcc > 0, ncc / torch.where(qcc > 0, qcc, torch.ones_like(qcc)), torch.zeros_like(qcc))
        au, ac = raw[:, 0], raw[:, 1]
        rates = torch.stack([au, ac, au / self.m_star, self.drops_per_embryo * au / self.m_star,
                             self.ac_n_factor * ac * nc_over_qc, -raw[:, 2], raw[:, 3]], 1)
        full = torch.zeros((qc.numel(), 7), dtype=torch.float64, device=qc.device)
        full[cells] = rates

        o0, o1, o2, o3, o4, o5, o6 = [v.reshape(shape) for v in full.unbind(1)]
        return o0, o1, o2, o3, o4, o5, o6, use_cloud.double().reshape(shape), use_rain.double().reshape(shape)


class PythonBackendEmulator:
    """What the python backend of components/emulators calls: infer(inputs, outputs)."""

    def __init__(self, model_file: str):
        self.model = SdmWarmRain(model_file)

    def infer(self, inputs: Dict[str, np.ndarray], outputs: Dict[str, np.ndarray]):
        with warnings.catch_warnings():
            # (the inputs are read-only arrays, and are only read)
            warnings.filterwarnings('ignore', message='The given NumPy array is not writable')
            args = [torch.from_numpy(inputs[n]) for n in INPUTS]
        with torch.no_grad():
            res = self.model(*args)
        for name, v in zip(OUTPUTS, res):
            outputs[name][...] = v.numpy()


def create_emulator(config):
    # Each MPI rank uses its own threads only
    torch.set_num_threads(int(config.get('num_threads', 1)))
    return PythonBackendEmulator(config['model_path'])
