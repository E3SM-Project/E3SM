"""Tests for the example emulator's packing and output splitting.

That packing is the part any emulator has to get right: laying several named
fields out side by side, splitting the answer back over the declared outputs,
and the widths that have to agree for either to be right.  Needs numpy and
nothing else.
"""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from e3sm_emulator.context import Context
from e3sm_emulator.example import ExampleEmulator


def emulator(inputs, outputs, weight=None, bias=None, **config) -> ExampleEmulator:
    """An ExampleEmulator with an optional weight matrix set directly."""
    built = ExampleEmulator(
        {"inputs": inputs, "outputs": outputs, **config}, Context()
    )
    if weight is not None:
        built.weight = np.asarray(weight, dtype=np.float64)
        built.bias = np.zeros(built.weight.shape[1]) if bias is None else bias
    return built


class TestExampleEmulator(unittest.TestCase):
    def test_fields_are_packed_side_by_side_and_split_back(self):
        # T is one value per column, q is three: the model sees them laid out
        # as [ncol, 4] in the declared order.  This W maps that to
        # [2*T, -q0, -q1], which the declared outputs consume as widths 1 and 2.
        W = np.zeros((4, 3))
        W[0, 0] = 2.0
        W[1, 1] = -1.0
        W[2, 2] = -1.0

        T = np.arange(4.0)
        q = np.arange(12.0).reshape(4, 3)
        dT = np.zeros(4)
        dq = np.zeros((4, 2))

        emulator(["T", "q"], ["dT", "dq"], weight=W).infer(
            {"T": T, "q": q}, {"dT": dT, "dq": dq}
        )

        np.testing.assert_allclose(dT, T * 2.0)
        np.testing.assert_allclose(dq, q[:, :2] * -1.0)

    def test_the_affine_default_writes_through_to_the_caller(self):
        # No weights: the point is that the result lands in the array E3SM
        # handed us, which is what the zero-copy view buys.
        T = np.arange(4.0)
        dT = np.zeros(4)
        emulator(["T"], ["dT"], scale=3.0, offset=1.0).infer({"T": T}, {"dT": dT})
        np.testing.assert_allclose(dT, T * 3.0 + 1.0)

    def test_weights_load_from_an_npz(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "weights.npz")
            np.savez(path, W=np.full((1, 1), 5.0), b=np.zeros(1))
            built = ExampleEmulator(
                {"inputs": ["T"], "outputs": ["dT"], "model_path": path}, Context()
            )
        T = np.arange(4.0)
        dT = np.zeros(4)
        built.infer({"T": T}, {"dT": dT})
        np.testing.assert_allclose(dT, T * 5.0)

    def test_a_rank_with_no_columns_does_nothing(self):
        # Normal on a large layout, and the reason this is not simply left to
        # numpy: reshape(0, -1) cannot infer the trailing extent, so every
        # packing step would raise.  A column-local model communicates
        # nothing, so there is no collective this rank fails to reach.
        emulator(["T"], ["dT"]).infer({"T": np.zeros(0)}, {"dT": np.zeros(0)})

    def test_a_model_that_returns_too_little_is_named(self):
        built = emulator(["T"], ["dT", "dq"], weight=np.ones((1, 1)))
        with self.assertRaises(ValueError) as caught:
            built.infer({"T": np.zeros(4)}, {"dT": np.zeros(4), "dq": np.zeros(4)})
        self.assertIn("'dq'", str(caught.exception))

    def test_a_model_that_returns_too_much_is_caught(self):
        # Silently dropping the tail would leave a field the namelist never
        # mentioned quietly unused.
        built = emulator(["T"], ["dT"], weight=np.ones((1, 3)))
        with self.assertRaises(ValueError) as caught:
            built.infer({"T": np.zeros(4)}, {"dT": np.zeros(4)})
        self.assertIn("consume 1", str(caught.exception))

    def test_inputs_that_disagree_about_the_columns_are_refused(self):
        built = emulator(["T", "q"], ["dT"])
        with self.assertRaises(ValueError) as caught:
            built.infer({"T": np.zeros(4), "q": np.zeros(3)}, {"dT": np.zeros(4)})
        self.assertIn("disagree about the column count", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
