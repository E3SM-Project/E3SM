from unittest.case import TestCase
from unittest.mock import Mock

import numpy as np
from cdms2.axis import TransientAxis
from cdms2.tvariable import TransientVariable

from e3sm_diags.driver.annual_cycle_zonal_mean_driver import _create_annual_cycle


class TestCreateAnnualCycle(TestCase):
    def test_returns_annual_cycle_for_a_dataset_variable(self):
        # Mock a Dataset object and get_climo_variable()
        dataset_mock = Mock()
        dataset_mock.get_climo_variable.return_value = TransientVariable(
            data=np.zeros((2, 2)),
            attributes={"id": "PRECNT", "long_name": "long_name", "units": "units"},
            axes=[
                TransientAxis(np.zeros(2), id="latitude"),
                TransientAxis(np.zeros(2), id="longitude"),
            ],
        )

        # Generate expected and result
        expected = TransientVariable(
            data=np.zeros((12, 2, 2)),
            attributes={"id": "PRECNT", "long_name": "long_name", "units": "units"},
            axes=[
                TransientAxis(
                    np.arange(1, 13),
                    id="time",
                    attributes={"axis": "T", "calendar": "gregorian"},
                ),
                TransientAxis(np.zeros(2), id="latitude"),
                TransientAxis(np.zeros(2), id="longitude"),
            ],
        )
        result = _create_annual_cycle(dataset_mock, variable="PRECNT")

        # Check data are equal
        np.array_equal(result.data, expected.data)

        # Check attributes are equal. Must delete "name" attribute since they differ.
        result.deleteattribute("name")
        expected.deleteattribute("name")
        self.assertDictEqual(result.attributes, expected.attributes)

        # Check time, latitude, and longitude axes are equal
        np.array_equal(result.getAxis(0)[:], expected.getAxis(0)[:])
        np.array_equal(result.getLatitude()[:], expected.getLatitude()[:])
        np.array_equal(result.getLongitude()[:], expected.getLongitude()[:])
