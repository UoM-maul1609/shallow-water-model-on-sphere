"""Utilities for periodic longitude grids.

The model now writes N distinct longitude samples over [0, 2*pi), rather than
including both 0 and 2*pi. These helpers also accept legacy files containing a
duplicated periodic endpoint.
"""

import numpy as np


def has_duplicated_endpoint(lons):
    """Return True when the first and last longitude are separated by ~2*pi."""
    lons = np.asarray(lons, dtype=float).reshape(-1)
    if lons.size < 2:
        return False
    diffs = np.diff(lons)
    step = np.median(np.abs(diffs)) if diffs.size else 0.0
    tol = max(1.0e-10, step * 1.0e-6)
    return np.isclose(abs(lons[-1] - lons[0]), 2.0 * np.pi, atol=tol, rtol=0.0)


def strip_duplicated_endpoint(lons, *fields):
    """Return native periodic samples, stripping the last longitude from legacy files.

    Field longitude is assumed to be the last axis.
    """
    lons = np.asarray(lons)
    duplicate = has_duplicated_endpoint(lons)
    lons_native = lons[:-1] if duplicate else lons
    fields_native = []
    for field in fields:
        arr = np.asanyarray(field)
        if arr.shape[-1] != lons.size:
            raise ValueError(
                f"longitude axis mismatch: field has {arr.shape[-1]} samples, "
                f"longitude has {lons.size}"
            )
        fields_native.append(arr[..., :-1] if duplicate else arr)
    return (lons_native, *fields_native)


def add_cyclic_longitude(lons, field):
    """Append a temporary cyclic longitude/data column for map plotting only."""
    lons_native, field_native = strip_duplicated_endpoint(lons, field)
    if lons_native.size < 2:
        return lons_native, field_native
    direction = np.sign(np.median(np.diff(lons_native)))
    if direction == 0:
        direction = 1.0
    lons_plot = np.concatenate(
        [lons_native, [lons_native[0] + direction * 2.0 * np.pi]]
    )
    field_plot = np.concatenate([field_native, field_native[..., 0:1]], axis=-1)
    return lons_plot, field_plot
