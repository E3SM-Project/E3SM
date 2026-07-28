"""Tests for the python emulator package.

Run directly (``python3 test_e3sm_emulator.py``) or under ctest.  Needs numpy
and nothing else: no MPI, no launcher, no checkpoint.
"""

from __future__ import annotations

import unittest

import numpy as np

from e3sm_emulator.bridge import create_emulator
from e3sm_emulator.context import Context


class TestContext(unittest.TestCase):
    def test_from_dict_matches_the_cxx_payload(self):
        # The keys here are exactly what python_inference_backend.cpp writes.
        context = Context.from_dict(
            {
                "rank": 2,
                "world_size": 4,
                "node_name": "nid001",
                "nx": 8,
                "ny": 6,
                "num_global_cols": 48,
                "col_gids": np.array([1, 5, 9], dtype=np.int32),
                "lat": np.array([-45.0, 0.0, 45.0]),
                "lon": np.array([0.0, 120.0, 240.0]),
            }
        )
        self.assertEqual(context.num_local_cols, 3)
        self.assertFalse(context.is_root)
        self.assertEqual(context.col_gids.dtype, np.int64)
        self.assertIn("3 of 48 columns", context.describe())

    def test_defaults_describe_a_serial_run(self):
        context = Context.from_dict({})
        self.assertTrue(context.is_root)
        self.assertEqual(context.world_size, 1)
        self.assertEqual(context.num_local_cols, 0)


class TestBridge(unittest.TestCase):
    def test_unknown_emulator_names_the_alternatives(self):
        with self.assertRaises(ValueError) as caught:
            create_emulator({"emulator": "nope"})
        message = str(caught.exception)
        self.assertIn("example", message)
        self.assertIn("python_module", message)

    def test_the_default_emulator_is_built(self):
        built = create_emulator({"context": {}, "inputs": ["T"], "outputs": ["dT"]})
        self.assertTrue(hasattr(built, "infer"))


if __name__ == "__main__":
    unittest.main()
