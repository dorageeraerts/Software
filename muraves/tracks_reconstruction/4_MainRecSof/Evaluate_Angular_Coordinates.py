from __future__ import annotations

import math


def TrackAngularCoordinates(slope_xy: float, slope_xz: float, X0: float, X2: float) -> list[float]:
    """Python port of TrackAngularCoordinates from EvaluateAngularCoordinates.cc."""
    dz = slope_xz * (X0 - X2)
    dx = X2 - X0
    dy = slope_xy * dx

    theta = math.atan(abs(dz) / math.sqrt(dx * dx + dy * dy)) * (180.0 / math.pi)

    if dz > 0:
        phi = math.atan(slope_xy) * (180.0 / math.pi) + 180.0
    else:
        if slope_xy > 0:
            phi = math.atan(slope_xy) * 180.0 / math.pi
        else:
            phi = 360.0 + math.atan(slope_xy) * 180.0 / math.pi

    return [theta, phi]


# Pythonic alias
track_angular_coordinates = TrackAngularCoordinates
