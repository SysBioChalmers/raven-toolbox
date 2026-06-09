"""The visualization subpackage is a documented stub (not yet implemented)."""
import pytest


def test_visualization_stub_raises_not_implemented():
    import raven_python.visualization as viz

    with pytest.raises(NotImplementedError, match="not implemented yet"):
        viz.draw_pathway()  # attribute access triggers __getattr__ and raises
