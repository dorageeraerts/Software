#!/usr/bin/env python3
"""Plot the directional geometrical acceptance of the ROSSO hodoscope.

The script evaluates the common overlap of equal rectangular sensitive planes.
For every direction (theta, phi), it calculates for each plane pair i,j

    Delta y_ij = Delta x_ij tan(phi)
    Delta z_ij = Delta x_ij tan(theta) / cos(phi)
    w_y,ij = [l_y - |Delta y_ij - O_y,ij|]_+
    w_z,ij = [l_z - |Delta z_ij - O_z,ij|]_+

and then uses the most restrictive pair in each projection:

    w_y = min(w_y,ij),  w_z = min(w_z,ij),  S_eff = w_y w_z.

No additional projected-area factor such as |cos(theta) cos(phi)| is applied.
With l_y = l_z = 1 m, S_eff in m^2 is numerically equal to g = S_eff/S.
"""

import argparse
import math
from dataclasses import dataclass
from itertools import combinations


# -----------------------------------------------------------------------------
# User-adjustable defaults
# -----------------------------------------------------------------------------

PHI_MIN_DEG = 0.0
PHI_MAX_DEG = 360.0
THETA_MIN_DEG = 0.0
THETA_MAX_DEG = 90.0
DPHI_DEG = 0.5
DTHETA_DEG = 0.5

SENSITIVE_LY_M = 1.0
SENSITIVE_LZ_M = 1.0

# Station-centre coordinates in metres. The origin is at the centre of plane 2.
# The z values follow from the measured base-height differences
# (300, 255, 211, and 0 mm for planes 1--4). Update these values here if a
# later geometry survey gives more precise station-centre coordinates.
PLANE_COORDINATES_M = {
    1: (-0.264, 0.0, +0.045),
    2: (+0.000, 0.0, +0.000),
    3: (+0.258, 0.0, -0.044),
    4: (+1.482, 0.0, -0.255),
}

PLANE_SETS = {
    "Four-plane acceptance (planes 1--4)": (1, 2, 3, 4),
    "Three-plane acceptance (planes 1--3)": (1, 2, 3),
}

PLOT_QUANTITY = "Seff"  # "Seff" for m^2 or "g" for the dimensionless factor
OUTPUT_PREFIX = "ROSSO_directional_geometrical_acceptance"


@dataclass(frozen=True)
class Plane:
    index: int
    x: float
    y: float
    z: float


PLANES = {
    index: Plane(index, *coordinates)
    for index, coordinates in PLANE_COORDINATES_M.items()
}


def positive_part(value):
    """Return [value]_+ = max(value, 0)."""

    return max(value, 0.0)


def pair_overlap_widths(theta_rad, phi_rad, plane_i, plane_j, ly, lz):
    """Return (w_y,ij, w_z,ij) for one pair of planes."""

    cos_phi = math.cos(phi_rad)
    if abs(cos_phi) < 1.0e-12:
        return 0.0, 0.0

    delta_x = plane_j.x - plane_i.x
    delta_y = delta_x * math.tan(phi_rad)
    delta_z = delta_x * math.tan(theta_rad) / cos_phi

    offset_y = plane_j.y - plane_i.y
    offset_z = plane_j.z - plane_i.z

    width_y = positive_part(ly - abs(delta_y - offset_y))
    width_z = positive_part(lz - abs(delta_z - offset_z))
    return width_y, width_z


def effective_area(theta_rad, phi_rad, plane_indices, ly=SENSITIVE_LY_M, lz=SENSITIVE_LZ_M):
    """Return (S_eff, w_y, w_z) for a selected set of tracking planes."""

    common_width_y = ly
    common_width_z = lz

    for index_i, index_j in combinations(plane_indices, 2):
        width_y, width_z = pair_overlap_widths(
            theta_rad, phi_rad, PLANES[index_i], PLANES[index_j], ly, lz
        )
        common_width_y = min(common_width_y, width_y)
        common_width_z = min(common_width_z, width_z)

        if common_width_y == 0.0 or common_width_z == 0.0:
            return 0.0, common_width_y, common_width_z

    return common_width_y * common_width_z, common_width_y, common_width_z


def number_of_bins(axis_min, axis_max, step):
    """Return an integer bin count and reject incompatible ranges."""

    if axis_max <= axis_min:
        raise ValueError("The axis maximum must be larger than the minimum.")
    if step <= 0.0:
        raise ValueError("The angular bin width must be positive.")

    unrounded = (axis_max - axis_min) / step
    rounded = int(round(unrounded))
    if not math.isclose(unrounded, rounded, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError("Each angular range must be an integer multiple of its bin width.")
    return rounded


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot three-plane and four-plane ROSSO geometrical-acceptance maps with PyROOT."
    )
    parser.add_argument("--phi-min", type=float, default=PHI_MIN_DEG)
    parser.add_argument("--phi-max", type=float, default=PHI_MAX_DEG)
    parser.add_argument("--theta-min", type=float, default=THETA_MIN_DEG)
    parser.add_argument("--theta-max", type=float, default=THETA_MAX_DEG)
    parser.add_argument("--dphi", type=float, default=DPHI_DEG)
    parser.add_argument("--dtheta", type=float, default=DTHETA_DEG)
    parser.add_argument("--quantity", choices=("Seff", "g"), default=PLOT_QUANTITY)
    parser.add_argument("--output", default=OUTPUT_PREFIX)
    return parser.parse_args()


def make_histogram(ROOT, name, plane_indices, args):
    n_phi = number_of_bins(args.phi_min, args.phi_max, args.dphi)
    n_theta = number_of_bins(args.theta_min, args.theta_max, args.dtheta)

    histogram = ROOT.TH2D(
        name,
        "",
        n_phi,
        args.phi_min,
        args.phi_max,
        n_theta,
        args.theta_min,
        args.theta_max,
    )
    histogram.SetDirectory(0)

    full_area = SENSITIVE_LY_M * SENSITIVE_LZ_M
    for bin_phi in range(1, n_phi + 1):
        phi_deg = histogram.GetXaxis().GetBinCenter(bin_phi)
        phi_rad = math.radians(phi_deg)

        for bin_theta in range(1, n_theta + 1):
            theta_deg = histogram.GetYaxis().GetBinCenter(bin_theta)
            theta_rad = math.radians(theta_deg)
            seff, _, _ = effective_area(theta_rad, phi_rad, plane_indices)
            value = seff if args.quantity == "Seff" else seff / full_area
            histogram.SetBinContent(bin_phi, bin_theta, value)

    z_title = "Geometrical Factor" if args.quantity == "Seff" else "g(#theta,#phi)"
    histogram.GetXaxis().SetTitle("Azimuth #phi [deg]")
    histogram.GetYaxis().SetTitle("Elevation #theta [deg]")
    histogram.GetZaxis().SetTitle(z_title)
    histogram.GetZaxis().SetTitleOffset(1.1)
    histogram.SetMinimum(0.0)
    histogram.SetMaximum(full_area if args.quantity == "Seff" else 1.0)
    histogram.SetContour(100)
    return histogram


def configure_style(ROOT):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetNumberContours(100)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelSize(0.042, "XYZ")
    ROOT.gStyle.SetTitleSize(0.048, "XYZ")


def main():
    args = parse_arguments()

    try:
        import ROOT
    except ImportError as error:
        raise SystemExit(
            "PyROOT is required to draw the maps. Run this script in an environment "
            "where `python3 -c 'import ROOT'` succeeds."
        ) from error

    ROOT.gROOT.SetBatch(True)
    configure_style(ROOT)

    histograms = []
    for map_index, (label, plane_indices) in enumerate(PLANE_SETS.items(), start=1):
        histogram = make_histogram(ROOT, f"h_acceptance_{map_index}", plane_indices, args)
        histograms.append((label, plane_indices, histogram))

    canvas = ROOT.TCanvas("c_acceptance", "Directional geometrical acceptance", 1200, 1050)
    canvas.Divide(1, 2, 0.0, 0.0)
    labels = []

    for pad_index, (label, plane_indices, histogram) in enumerate(histograms, start=1):
        pad = canvas.cd(pad_index)
        pad.SetLeftMargin(0.11)
        pad.SetRightMargin(0.15)
        pad.SetTopMargin(0.08)
        pad.SetBottomMargin(0.12)
        pad.SetTicks(1, 1)
        histogram.Draw("COLZ")

        text = ROOT.TLatex()
        text.SetNDC(True)
        text.SetTextFont(42)
        text.SetTextSize(0.050)
        text.DrawLatex(0.12, 0.94, label)
        labels.append(text)

    canvas.SaveAs(f"{args.output}.pdf")


    print(f"Saved {args.output}.pdf")


if __name__ == "__main__":
    main()
