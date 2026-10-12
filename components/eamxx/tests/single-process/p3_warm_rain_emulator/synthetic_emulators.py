"""Synthetic warm-rain emulators of test p3_warm_rain_emulator.

They follow the contract of P3's warm-rain emulator (eamxx_p3_warm_rain_emulator.cpp),
and need numpy only. warm_rain_emulator_model_path names the one to use:
  nowhere: never applies (masks are 0), with rates that would change P3 if they
           were used. P3 must then be BFB with P3 without it.
  no_rain: applies everywhere, with zero rates: the warm-rain collision processes
           are off.
  drain:   applies everywhere, and converts all the cloud water it is given to
           rain in less than a time step (of at least 50 s).
"""
import numpy as np

RATES = ('qc2qr_autoconv_tend', 'qc2qr_accret_tend', 'ncautr', 'nc2nr_autoconv_tend',
         'nc_accret_tend', 'nc_selfcollect_tend', 'nr_selfcollect_tend')
MASKS = ('use_cloud', 'use_rain')


class Synthetic:
    def __init__(self, kind):
        if kind not in ('nowhere', 'no_rain', 'drain'):
            raise ValueError(f"synthetic_emulators: unknown model '{kind}'")
        self.kind = kind

    def infer(self, inputs, outputs):
        nowhere = self.kind == 'nowhere'
        for name in RATES:
            outputs[name][...] = 1.0e3 if nowhere else 0.0
        if nowhere:
            outputs['nc_selfcollect_tend'][...] = -1.0e3
        if self.kind == 'drain':
            outputs['qc2qr_autoconv_tend'][...] = np.maximum(inputs['qc'], 0.0) / 50.0
        for name in MASKS:
            outputs[name][...] = 0.0 if nowhere else 1.0


def create_emulator(config):
    return Synthetic(config['model_path'])
