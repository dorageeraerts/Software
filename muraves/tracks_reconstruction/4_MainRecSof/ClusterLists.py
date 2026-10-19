from __future__ import annotations

from dataclasses import dataclass, field
import random
from typing import Sequence

from reco_config import get_reco_config


@dataclass
class ClusterCollection:
    ClustersEnergy: list[float] = field(default_factory=list)
    ClustersPositions: list[float] = field(default_factory=list)
    StripsEnergy: list[list[float]] = field(default_factory=list)
    StripsPositions: list[list[float]] = field(default_factory=list)
    StripsID: list[list[float]] = field(default_factory=list)
    TimeExpansions: list[float] = field(default_factory=list)
    ClustersSize: list[int] = field(default_factory=list)


class DeterministicSmearingRNG:
    """Cross-language RNG for reproducible single-strip smearing."""

    def __init__(self, seed: int):
        state = int(seed) & 0xFFFFFFFF
        self._state = state if state != 0 else 0x6D2B79F5

    def _next_u32(self) -> int:
        state = self._state
        state ^= (state << 13) & 0xFFFFFFFF
        state ^= (state >> 17) & 0xFFFFFFFF
        state ^= (state << 5) & 0xFFFFFFFF
        self._state = state & 0xFFFFFFFF
        return self._state

    def randint(self, low: int, high: int) -> int:
        if high < low:
            raise ValueError("high must be >= low")
        return low + (self._next_u32() % (high - low + 1))


def ClusterPosition(stripDeposits: Sequence[float], stripPos: Sequence[float]) -> float:
    """Return weighted cluster position from strip deposits and positions."""
    numerator = 0.0
    denominator = 0.0
    for deposit, position in zip(stripDeposits, stripPos):
        numerator += deposit * position
        denominator += deposit
    return numerator / denominator


def ClusterEnergy(stripDeposits: Sequence[float]) -> float:
    """Return total cluster energy as sum of strip deposits."""
    return float(sum(stripDeposits))


def SortIndices(referenceVector: Sequence[float]) -> list[int]:
    """Sort indices by decreasing values of the reference vector."""
    return sorted(range(len(referenceVector)), key=lambda idx: referenceVector[idx], reverse=True)


def CreateClusterList(
    Deposits: Sequence[float],
    EnergyThreshold_clusterStrip: float,
    EnergyThreshold_singleStrip: float,
    AdStripsThEnergy_singleStripCl: float,
    Texp1: float,
    Texp2: float,
    TriggerMask1: Sequence[float],
    TriggerMask2: Sequence[float],
    smearing_rng: DeterministicSmearingRNG | None = None,
    cluster_cfg: dict | None = None,
) -> ClusterCollection:
    """
    Python port of CreateClusterList from ClusterLists.cc.

    The implementation intentionally mirrors the C++ logic, including cluster
    ordering by energy and the single-strip random position smearing.
    """

    cfg = cluster_cfg if cluster_cfg is not None else get_reco_config()["cluster_lists"]
    max_strip_energy = float(cfg["max_strip_energy"])
    first_strip_pos = float(cfg["first_strip_pos"])
    adjacent_strips_distance = float(cfg["adjacent_strips_distance"])
    n_strips = int(cfg["n_strips"])
    min_energy_trigger_mask = float(cfg["min_energy_trigger_mask"])
    split_board_index = int(cfg["split_board_index"])
    smearing_scale = int(cfg["single_strip_smearing_scale"])

    cluster_strips_deposits: list[float] = []
    cluster_strips_pos: list[float] = []
    cluster_strips_id: list[float] = []

    cluster_positions: list[float] = []
    cluster_energies: list[float] = []
    cluster_sizes: list[int] = []

    clusters_strips_deposits: list[list[float]] = []
    clusters_strips_positions: list[list[float]] = []
    clusters_strips_id: list[list[float]] = []

    strip_index: list[int] = []
    clu_en = 0.0
    clu_pos = 0.0
    clu_size = 0

    # Preserve the old behavior unless a deterministic cross-language RNG is provided.
    nondeterministic_rng = random.SystemRandom() if smearing_rng is None else None

    for st in range(n_strips):
        deposit = float(Deposits[st])

        if deposit < EnergyThreshold_clusterStrip:
            cluster_is_on = False
        else:
            cluster_is_on = True

        if deposit > max_strip_energy:
            cluster_is_on = False

        if deposit > min_energy_trigger_mask:
            # condition regarding strips of the first borad (0-31) 
            if st < split_board_index:
                # if the strips is not in the trigger mask of the first board, there is no cluster.
                if st not in TriggerMask1:
                    cluster_is_on = False
            # condition regarding strips belonging to the second board (32-63)
            else:
                # if the strips is not in the trigger mask of the second board, there is no cluster.
                if (st - split_board_index) not in TriggerMask2:
                    cluster_is_on = False
        
        # (?) In the C++ code, there is an additional condition that checks if the strip index is the last one (n_strips - 1) and if so, it turns off the cluster.
        if st == (n_strips - 1):
            cluster_is_on = False

        if cluster_is_on:
            strip_position = st * adjacent_strips_distance
            cluster_strips_deposits.append(deposit)
            cluster_strips_pos.append(first_strip_pos + strip_position)
            cluster_strips_id.append(float(st))
            clu_en += deposit
            # cluster position is calculated as a weighted average of the strip positions, where the weights are the strip deposits.
            clu_pos += (strip_position + first_strip_pos) * deposit 
            clu_size += 1
            strip_index.append(st)
        else:
            # processing the case of single-strip cluster, which is a special case in the C++ code.
            if clu_size == 1:
                if EnergyThreshold_singleStrip < clu_en < max_strip_energy:
                    base_strip = strip_index[0]
                    # Checking strips adjacent to the single-strip cluster to see if they can be included in the cluster based on their energy deposits.
                    # Strip adjacent to the right (base_strip + 1)
                    if (
                        base_strip < (n_strips - 1)
                        and Deposits[base_strip + 1] > AdStripsThEnergy_singleStripCl
                        and Deposits[base_strip + 1] < max_strip_energy
                    ):
                        post_strip_position = (base_strip + 1) * adjacent_strips_distance
                        clu_pos += Deposits[base_strip + 1] * (first_strip_pos + post_strip_position)
                        clu_en += Deposits[base_strip + 1]
                        clu_size += 1
                        cluster_strips_deposits.append(float(Deposits[base_strip + 1]))
                        cluster_strips_pos.append(first_strip_pos + post_strip_position)
                        cluster_strips_id.append(float(base_strip + 1))
                    # Strip adjacent to the left (base_strip - 1)
                    if (
                        base_strip > 0
                        and Deposits[base_strip - 1] > AdStripsThEnergy_singleStripCl
                        and Deposits[base_strip - 1] < max_strip_energy
                    ):
                        pre_strip_position = (base_strip - 1) * adjacent_strips_distance
                        clu_pos += Deposits[base_strip - 1] * (first_strip_pos + pre_strip_position)
                        clu_en += Deposits[base_strip - 1]
                        cluster_strips_deposits.append(float(Deposits[base_strip - 1]))
                        cluster_strips_pos.append(first_strip_pos + pre_strip_position)
                        cluster_strips_id.append(float(base_strip - 1))
                        clu_size += 1
                else:
                    clu_pos = 0.0
                    clu_en = 0.0
                    clu_size = 0

            if clu_size > 0:
                # The cluster position is calculated as the weighted average of the strip positions, where the weights are the strip deposits.
                # Here, we divide the accumulated position by the total energy to get the final cluster position to complete the weighted average calculation.
                cluster_position = clu_pos / clu_en
                if clu_size == 1:
                    # Keep C++ granularity: integer millimeter steps after 1/1000 scaling.
                    range_from = int(-smearing_scale * adjacent_strips_distance / 2)
                    range_to = int(smearing_scale * adjacent_strips_distance / 2)
                    rng = smearing_rng if smearing_rng is not None else nondeterministic_rng
                    cluster_position += rng.randint(range_from, range_to) / float(smearing_scale)

                cluster_positions.append(cluster_position)
                cluster_energies.append(clu_en)
                cluster_sizes.append(clu_size)

                ordered_strip_indices = SortIndices(cluster_strips_deposits)
                sorted_strip_ids = [cluster_strips_id[idx] for idx in ordered_strip_indices]
                sorted_strip_positions = [cluster_strips_pos[idx] for idx in ordered_strip_indices]
                sorted_strip_deposits = [cluster_strips_deposits[idx] for idx in ordered_strip_indices]

                clusters_strips_deposits.append(sorted_strip_deposits)
                clusters_strips_positions.append(sorted_strip_positions)
                clusters_strips_id.append(sorted_strip_ids)

            strip_index.clear()
            cluster_strips_deposits.clear()
            cluster_strips_pos.clear()
            cluster_strips_id.clear()
            clu_en = 0.0
            clu_pos = 0.0
            clu_size = 0

    ordered_cluster_indices = SortIndices(cluster_energies)

    sorted_cluster_positions: list[float] = []
    sorted_cluster_energies: list[float] = []
    sorted_cluster_size: list[int] = []
    sorted_clusters_strips_positions: list[list[float]] = []
    sorted_clusters_strips_deposits: list[list[float]] = []
    sorted_clusters_strips_id: list[list[float]] = []
    texp_cluster: list[float] = []

    split_boundary_low = first_strip_pos + (split_board_index - 1) * adjacent_strips_distance
    split_boundary_high = first_strip_pos + split_board_index * adjacent_strips_distance

    for idx in ordered_cluster_indices:
        pos = cluster_positions[idx]
        sorted_cluster_positions.append(pos)
        sorted_cluster_energies.append(cluster_energies[idx])
        sorted_cluster_size.append(cluster_sizes[idx])
        sorted_clusters_strips_positions.append(clusters_strips_positions[idx])
        sorted_clusters_strips_deposits.append(clusters_strips_deposits[idx])
        sorted_clusters_strips_id.append(clusters_strips_id[idx])

        if pos <= split_boundary_low:
            texp = Texp1
        else:
            texp = Texp2

        if split_boundary_low < pos < split_boundary_high:
            texp = Texp1 if Texp1 > Texp2 else Texp2

        texp_cluster.append(texp)

    return ClusterCollection(
        ClustersEnergy=sorted_cluster_energies,
        ClustersPositions=sorted_cluster_positions,
        StripsEnergy=sorted_clusters_strips_deposits,
        StripsPositions=sorted_clusters_strips_positions,
        StripsID=sorted_clusters_strips_id,
        TimeExpansions=texp_cluster,
        ClustersSize=sorted_cluster_size,
    )


# Pythonic aliases
cluster_position = ClusterPosition
cluster_energy = ClusterEnergy
sort_indices = SortIndices
create_cluster_list = CreateClusterList
