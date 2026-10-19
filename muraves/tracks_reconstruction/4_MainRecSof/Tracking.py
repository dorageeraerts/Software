from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Sequence

from reco_config import get_reco_config


@dataclass
class TracksCollection:
    intercepts_3p: list[float] = field(default_factory=list)
    slopes_3p: list[float] = field(default_factory=list)
    chiSquares_3p: list[float] = field(default_factory=list)
    cluster_index_1: list[int] = field(default_factory=list)
    cluster_index_2: list[int] = field(default_factory=list)
    cluster_index_3: list[int] = field(default_factory=list)
    Ntracks_4p: list[int] = field(default_factory=list)
    position_c1: list[float] = field(default_factory=list)
    position_c2: list[float] = field(default_factory=list)
    position_c3: list[float] = field(default_factory=list)
    residue_c1: list[float] = field(default_factory=list)
    residue_c2: list[float] = field(default_factory=list)
    residue_c3: list[float] = field(default_factory=list)
    TrackEnergy_3p: list[float] = field(default_factory=list)
    IsInTrack_clusters1: list[int] = field(default_factory=list)
    IsInTrack_clusters2: list[int] = field(default_factory=list)
    IsInTrack_clusters3: list[int] = field(default_factory=list)
    IsInTrack_clusters1_4p: list[int] = field(default_factory=list)
    IsInTrack_clusters2_4p: list[int] = field(default_factory=list)
    IsInTrack_clusters3_4p: list[int] = field(default_factory=list)
    Plane4th_isIntercepted: list[int] = field(default_factory=list)
    positions_c4: list[list[float]] = field(default_factory=list)
    intercept_4p: list[list[float]] = field(default_factory=list)
    slope_4p: list[list[float]] = field(default_factory=list)
    chiSquares_4p: list[list[float]] = field(default_factory=list)
    displacement_p4: list[list[float]] = field(default_factory=list)
    cluster_indices_4: list[list[float]] = field(default_factory=list)
    residue_c1_p4: list[list[float]] = field(default_factory=list)
    residue_c2_p4: list[list[float]] = field(default_factory=list)
    residue_c3_p4: list[list[float]] = field(default_factory=list)
    residue_c4_p4: list[list[float]] = field(default_factory=list)
    TrackEnergy_4p: list[list[float]] = field(default_factory=list)
    ScatteringAngles: list[list[float]] = field(default_factory=list)
    ExpectedPosition_OnPlane4th: list[float] = field(default_factory=list)
    Track_3p_to_4p_index: list[int] = field(default_factory=list)
    Track_3p_ExpectedRes_p2: list[float] = field(default_factory=list)
    BestEnergy: float = 0.0
    BestChi: float = 10000.0
    BestEnergy_index: int = -1
    BestChi_index: int = -1


def _linear_fit_with_known_sigma(x: Sequence[float], y: Sequence[float], sigma: float) -> tuple[float, float, float, float, float]:
    """
    Weighted linear fit for y = intercept + slope * x with uniform y uncertainty.

    Returns (intercept, slope, chi_square, intercept_error, slope_error).
    """
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    x_mean = sum(x) / n
    y_mean = sum(y) / n

    sxx = sum((xi - x_mean) * (xi - x_mean) for xi in x)
    sxy = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))

    if sxx == 0.0:
        slope = 0.0
    else:
        slope = sxy / sxx

    intercept = y_mean - slope * x_mean

    if sigma > 0:
        chi_square = sum(((yi - (intercept + slope * xi)) / sigma) ** 2 for xi, yi in zip(x, y))
    else:
        chi_square = sum((yi - (intercept + slope * xi)) ** 2 for xi, yi in zip(x, y))

    if sxx > 0.0 and sigma > 0:
        slope_error = sigma / math.sqrt(sxx)
        intercept_error = sigma * math.sqrt((1.0 / n) + (x_mean * x_mean / sxx))
    else:
        slope_error = 0.0
        intercept_error = 0.0

    return intercept, slope, chi_square, intercept_error, slope_error


def MakeTracks(
    clusters_1: Sequence[float],
    clusters_2: Sequence[float],
    clusters_3: Sequence[float],
    clusters_4: Sequence[float],
    clustersEn_1: Sequence[float],
    clustersEn_2: Sequence[float],
    clustersEn_3: Sequence[float],
    clustersEn_4: Sequence[float],
    Texp_cl1: Sequence[float],
    Texp_cl2: Sequence[float],
    Texp_cl3: Sequence[float],
    Texp_cl4: Sequence[float],
    proximity_cut: float,
    X_pos: Sequence[float],
    Z_add: Sequence[float],
    sigma: float,
    tracking_cfg: dict | None = None,
) -> TracksCollection:
    """Python port of MakeTracks from Tracking.cc."""
    del Texp_cl4  # kept for signature compatibility with C++

    cfg = tracking_cfg if tracking_cfg is not None else get_reco_config()["tracking"]
    cluster_cfg = get_reco_config()["cluster_lists"]
    first_strip_pos = float(cfg["first_strip_pos"])
    adjacent_strips_distance = float(cfg["adjacent_strips_distance"])
    scattering_pi = float(cfg["scattering_pi"])
    n_strips = int(cluster_cfg["n_strips"])
    min_p4 = first_strip_pos
    max_p4 = first_strip_pos + (n_strips - 1) * adjacent_strips_distance

    intercept_3ptracks: list[float] = []
    slope_3ptracks: list[float] = []
    chiSquare_3ptracks: list[float] = []
    track_energy_3p: list[float] = []

    indices_clu1: list[int] = []
    indices_clu2: list[int] = []
    indices_clu3: list[int] = []

    position_c1: list[float] = []
    position_c2: list[float] = []
    position_c3: list[float] = []
    residue_c1: list[float] = []
    residue_c2: list[float] = []
    residue_c3: list[float] = []

    ntracks_4p: list[int] = []
    position_c4: list[list[float]] = []
    intercept_4ptracks: list[list[float]] = []
    slope_4ptracks: list[list[float]] = []
    chiSquare_4ptracks: list[list[float]] = []
    displacement_p4: list[list[float]] = []
    multiple_indices_clu4: list[list[float]] = []
    residue_c1_p4: list[list[float]] = []
    residue_c2_p4: list[list[float]] = []
    residue_c3_p4: list[list[float]] = []
    residue_c4_p4: list[list[float]] = []
    track_energies_4p: list[list[float]] = []
    scattering_angles: list[list[float]] = []

    expected_position_on_plane4th: list[float] = []
    plane4th_is_intercepted: list[int] = []
    tracks_3p_index: list[int] = []
    expected_position_res_p2: list[float] = []

    vect_of_vect_cl2: list[list[int]] = []
    vect_of_vect_cl3: list[list[int]] = []
    vect_of_vect_cl2_4p: list[list[int]] = []
    vect_of_vect_cl3_4p: list[list[int]] = []

    is_in_track_clusters1: list[int] = []
    is_in_track_clusters2: list[int] = []
    is_in_track_clusters3: list[int] = []
    is_in_track_clusters1_4p: list[int] = []
    is_in_track_clusters2_4p: list[int] = []
    is_in_track_clusters3_4p: list[int] = []

    best_chi_index = -1
    best_energy_index = -1
    best_chi = float("inf")
    best_energy = 0.0
    ntracks = 0

    for i in range(len(clusters_1)):
        ntrack_cl1 = 0
        ntrack_cl1_4p = 0

        p1 = clusters_1[i] + Z_add[0]
        x1 = X_pos[0]

        n_track_cl3_vector: list[int] = []
        n_track_cl3_vector_4p: list[int] = []

        for k in range(len(clusters_3)):
            ntrack_cl3 = 0
            ntrack_cl3_4p = 0

            p3 = clusters_3[k] + Z_add[2]
            x3 = X_pos[2]

            n_track_cl2_vector: list[int] = []
            n_track_cl2_vector_4p: list[int] = []

            for j in range(len(clusters_2)):
                x2 = X_pos[1]
                exp_m = (p1 - p3) / (x1 - x3)
                exp_q = p1 - ((p1 - p3) / (x1 - x3)) * x1
                exp_p2 = exp_m * x2 + exp_q

                p2 = clusters_2[j] + Z_add[1]
                exp_res_2 = p2 - exp_p2

                if abs(exp_res_2) < proximity_cut:
                    expected_position_res_p2.append(exp_res_2)
                    ntracks += 1
                    ntrack_cl3 += 1
                    ntrack_cl1 += 1

                    n_track_cl2_vector.append(1)
                    indices_clu1.append(i)
                    indices_clu2.append(j)
                    indices_clu3.append(k)

                    x_3p = [x1, x2, x3]
                    p_3p = [p1, p2, p3]
                    fit3p_intercept, fit3p_slope, fit3p_chi_square, intercept_error, slope_error = _linear_fit_with_known_sigma(
                        x_3p, p_3p, sigma
                    )

                    intercept_3ptracks.append(fit3p_intercept)
                    slope_3ptracks.append(fit3p_slope)
                    chiSquare_3ptracks.append(fit3p_chi_square)

                    position_c1.append(p1)
                    position_c2.append(p2)
                    position_c3.append(p3)
                    residue_c1.append(p1 - fit3p_intercept - fit3p_slope * x1)
                    residue_c2.append(p2 - fit3p_intercept - fit3p_slope * x2)
                    residue_c3.append(p3 - fit3p_intercept - fit3p_slope * x3)

                    track_energy_3p_single = clustersEn_1[i] + clustersEn_2[j] + clustersEn_3[k]
                    track_energy_3p.append(track_energy_3p_single)

                    if Texp_cl1[i] > 0 and Texp_cl2[j] and Texp_cl3[k]:
                        if best_chi > fit3p_chi_square:
                            best_chi = fit3p_chi_square
                            best_chi_index = ntracks - 1
                        if best_energy < track_energy_3p_single:
                            best_energy = track_energy_3p_single
                            best_energy_index = ntracks - 1

                    single_position_c4: list[float] = []
                    single_intercept_4ptracks: list[float] = []
                    single_slope_4ptracks: list[float] = []
                    single_chiSquare_4ptracks: list[float] = []
                    single_displacement_p4: list[float] = []
                    single_indices_clu4: list[float] = []
                    single_residue_c1_p4: list[float] = []
                    single_residue_c2_p4: list[float] = []
                    single_residue_c3_p4: list[float] = []
                    single_residue_c4_p4: list[float] = []
                    single_track_energy_4p: list[float] = []
                    theta_scatter: list[float] = []

                    x4 = X_pos[3]
                    exp_p4 = fit3p_intercept + x4 * fit3p_slope
                    _exp_p4_err = math.sqrt(intercept_error * intercept_error + x4 * x4 * slope_error * slope_error)

                    if (exp_p4 > (min_p4 - sigma)) and (exp_p4 < (max_p4 + sigma)):
                        expected_position_on_plane4th.append(exp_p4)
                        plane4th_is_intercepted.append(1)

                        if len(clusters_4) > 0:
                            ntrack_cl3_4p += 1
                            ntrack_cl1_4p += 1
                            n_track_cl2_vector_4p.append(1)
                            tracks_3p_index.append(ntracks - 1)
                        else:
                            n_track_cl2_vector_4p.append(0)

                        res_p4_allowed = [clusters_4[h] - exp_p4 for h in range(len(clusters_4))]
                        ntracks_4p.append(len(res_p4_allowed))

                        sorted_indices = sorted(
                            range(len(res_p4_allowed)),
                            key=lambda idx: abs(res_p4_allowed[idx]),
                        )

                        sorted_res_p4_allowed = [clusters_4[idx] for idx in sorted_indices]
                        single_indices_clu4 = [float(idx) for idx in sorted_indices]
                        single_displacement_p4 = [res_p4_allowed[idx] for idx in sorted_indices]

                        for cl4, p4 in enumerate(sorted_res_p4_allowed):
                            x_4p = [x1, x2, x3, x4]
                            p_4p = [p1, p2, p3, p4]
                            fit4p_intercept, fit4p_slope, fit4p_chi_square, _, _ = _linear_fit_with_known_sigma(
                                x_4p, p_4p, sigma
                            )

                            single_position_c4.append(p4)
                            single_intercept_4ptracks.append(fit4p_intercept)
                            single_slope_4ptracks.append(fit4p_slope)
                            single_chiSquare_4ptracks.append(fit4p_chi_square)
                            single_residue_c1_p4.append(p1 - fit4p_intercept - fit4p_slope * x1)
                            single_residue_c2_p4.append(p2 - fit4p_intercept - fit4p_slope * x2)
                            single_residue_c3_p4.append(p3 - fit4p_intercept - fit4p_slope * x3)
                            single_residue_c4_p4.append(p4 - fit4p_intercept - fit4p_slope * x4)

                            idx_for_energy = sorted_indices[cl4]
                            single_track_energy_4p.append(
                                clustersEn_1[i] + clustersEn_2[j] + clustersEn_3[k] + clustersEn_4[idx_for_energy]
                            )

                            side_a = math.sqrt(
                                (fit3p_intercept + x3 * fit3p_slope - fit3p_intercept - x4 * fit3p_slope)
                                * (fit3p_intercept + x3 * fit3p_slope - fit3p_intercept - x4 * fit3p_slope)
                                + (x3 - x4) * (x3 - x4)
                            )
                            side_b = math.sqrt(
                                (fit3p_intercept + x3 * fit3p_slope - p4)
                                * (fit3p_intercept + x3 * fit3p_slope - p4)
                                + (x3 - x4) * (x3 - x4)
                            )

                            denom = 2.0 * side_a * side_b
                            if denom == 0.0:
                                theta_scatter.append(0.0)
                            else:
                                cos_scatt = (
                                    (side_a * side_a)
                                    + (side_b * side_b)
                                    - (fit3p_intercept + x4 * fit3p_slope - p4)
                                    * (fit3p_intercept + x4 * fit3p_slope - p4)
                                ) / denom
                                cos_scatt = max(-1.0, min(1.0, cos_scatt))

                                displacement = single_displacement_p4[cl4]
                                if displacement == 0.0:
                                    sign = 0.0
                                else:
                                    sign = displacement / abs(displacement)
                                theta_scatter.append(sign * math.acos(cos_scatt) * (180.0 / scattering_pi))

                        scattering_angles.append(theta_scatter)
                        position_c4.append(single_position_c4)
                        intercept_4ptracks.append(single_intercept_4ptracks)
                        slope_4ptracks.append(single_slope_4ptracks)
                        chiSquare_4ptracks.append(single_chiSquare_4ptracks)
                        displacement_p4.append(single_displacement_p4)
                        multiple_indices_clu4.append(single_indices_clu4)
                        residue_c1_p4.append(single_residue_c1_p4)
                        residue_c2_p4.append(single_residue_c2_p4)
                        residue_c3_p4.append(single_residue_c3_p4)
                        residue_c4_p4.append(single_residue_c4_p4)
                        track_energies_4p.append(single_track_energy_4p)
                    else:
                        plane4th_is_intercepted.append(0)
                        ntracks_4p.append(0)
                        n_track_cl2_vector_4p.append(0)
                else:
                    n_track_cl2_vector.append(0)
                    n_track_cl2_vector_4p.append(0)

            n_track_cl3_vector.append(1 if ntrack_cl3 > 0 else 0)
            n_track_cl3_vector_4p.append(1 if ntrack_cl3_4p > 0 else 0)
            vect_of_vect_cl2.append(n_track_cl2_vector)
            vect_of_vect_cl2_4p.append(n_track_cl2_vector_4p)

        vect_of_vect_cl3.append(n_track_cl3_vector)
        vect_of_vect_cl3_4p.append(n_track_cl3_vector_4p)
        is_in_track_clusters1_4p.append(ntrack_cl1_4p)
        is_in_track_clusters1.append(ntrack_cl1)

    if len(vect_of_vect_cl3) > 0:
        for j in range(len(vect_of_vect_cl3[0])):
            n_t_cl = 0
            for i in range(len(vect_of_vect_cl3)):
                n_t_cl += vect_of_vect_cl3[i][j]
            is_in_track_clusters3.append(n_t_cl)

    if len(vect_of_vect_cl2) > 0:
        for j in range(len(vect_of_vect_cl2[0])):
            n_t_cl = 0
            for i in range(len(vect_of_vect_cl2)):
                n_t_cl += vect_of_vect_cl2[i][j]
            is_in_track_clusters2.append(n_t_cl)

    if len(vect_of_vect_cl3_4p) > 0:
        for j in range(len(vect_of_vect_cl3_4p[0])):
            n_t_cl = 0
            for i in range(len(vect_of_vect_cl3_4p)):
                n_t_cl += vect_of_vect_cl3_4p[i][j]
            is_in_track_clusters3_4p.append(n_t_cl)

    if len(vect_of_vect_cl2_4p) > 0:
        for j in range(len(vect_of_vect_cl2_4p[0])):
            n_t_cl = 0
            for i in range(len(vect_of_vect_cl2_4p)):
                n_t_cl += vect_of_vect_cl2_4p[i][j]
            is_in_track_clusters2_4p.append(n_t_cl)

    return TracksCollection(
        intercepts_3p=intercept_3ptracks,
        slopes_3p=slope_3ptracks,
        chiSquares_3p=chiSquare_3ptracks,
        cluster_index_1=indices_clu1,
        cluster_index_2=indices_clu2,
        cluster_index_3=indices_clu3,
        Ntracks_4p=ntracks_4p,
        position_c1=position_c1,
        position_c2=position_c2,
        position_c3=position_c3,
        residue_c1=residue_c1,
        residue_c2=residue_c2,
        residue_c3=residue_c3,
        TrackEnergy_3p=track_energy_3p,
        IsInTrack_clusters1=is_in_track_clusters1,
        IsInTrack_clusters2=is_in_track_clusters2,
        IsInTrack_clusters3=is_in_track_clusters3,
        IsInTrack_clusters1_4p=is_in_track_clusters1_4p,
        IsInTrack_clusters2_4p=is_in_track_clusters2_4p,
        IsInTrack_clusters3_4p=is_in_track_clusters3_4p,
        Plane4th_isIntercepted=plane4th_is_intercepted,
        positions_c4=position_c4,
        intercept_4p=intercept_4ptracks,
        slope_4p=slope_4ptracks,
        chiSquares_4p=chiSquare_4ptracks,
        displacement_p4=displacement_p4,
        cluster_indices_4=multiple_indices_clu4,
        residue_c1_p4=residue_c1_p4,
        residue_c2_p4=residue_c2_p4,
        residue_c3_p4=residue_c3_p4,
        residue_c4_p4=residue_c4_p4,
        TrackEnergy_4p=track_energies_4p,
        ScatteringAngles=scattering_angles,
        ExpectedPosition_OnPlane4th=expected_position_on_plane4th,
        Track_3p_to_4p_index=tracks_3p_index,
        Track_3p_ExpectedRes_p2=expected_position_res_p2,
        BestEnergy=best_energy,
        BestChi=best_chi,
        BestEnergy_index=best_energy_index,
        BestChi_index=best_chi_index,
    )


# Pythonic alias
make_tracks = MakeTracks
