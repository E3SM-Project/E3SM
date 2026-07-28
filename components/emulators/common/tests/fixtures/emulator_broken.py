"""An emulator whose factory fails, so the error path can be tested."""


def create_emulator(config):
    raise RuntimeError("deliberate failure from the test fixture")
