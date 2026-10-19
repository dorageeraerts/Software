from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from reco_config import get_reco_config


@dataclass
class EventStruct:
    boards: list[list[float]] = field(default_factory=list)
    TrMask_channels: list[list[float]] = field(default_factory=list)
    TrMask_strips: list[list[float]] = field(default_factory=list)
    TrMask_size: list[int] = field(default_factory=list)
    timeExp: list[float] = field(default_factory=list)
    timeStamp: float = 0.0


def _safe_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def ReadEvent(ADCline: str, stripIndices: Sequence[int], events_cfg: dict | None = None,) -> EventStruct:
    """Python port of ReadEvent from ReadEvent.cc."""
    adcdata_splitted = ADCline.rstrip("\n").split("\t")

    cfg = events_cfg if events_cfg is not None else get_reco_config()["read_event"]
    n_info_board = int(cfg["n_info_board"])
    n_channels = int(cfg["n_channels"])
    n_boards = int(cfg["n_boards"])
    time_stamp_index = int(cfg["time_stamp_index"])
    time_exp_offset = int(cfg["time_exp_offset"])
    trigger_mask_offset = int(cfg["trigger_mask_offset"])
    adc_offset = int(cfg["adc_offset"])
    trigger_mask_separator = str(cfg["trigger_mask_separator"])

    boards: list[list[float]] = []
    all_time_exp: list[float] = []

    boards_trigger_masks_channels: list[list[float]] = []
    boards_trigger_masks_strips: list[list[float]] = []
    trigger_mask_sizes: list[int] = []

    time_stamp = (
        _safe_float(adcdata_splitted[time_stamp_index])
        if len(adcdata_splitted) > time_stamp_index
        else 0.0
    )

    for n in range(n_boards):
        base_idx = n * n_info_board
        time_exp_idx = base_idx + time_exp_offset
        trmask_idx = base_idx + trigger_mask_offset
        
        # When it is read to be zero, it means that the board was not involved in the trigger, so we can safely set it to 0.0 in that case.
        # This variables is used later to apply a time cut on the event, namely to check if the time of the event is within a certain range from the trigger time. If the board was not involved in the trigger, it should not contribute to the timing of the event, so setting it to 0.0 is a reasonable default that allows us to ignore it in timing calculations.
        time_exp = _safe_float(adcdata_splitted[time_exp_idx]) if len(adcdata_splitted) > time_exp_idx else 0.0
        all_time_exp.append(time_exp)

        trmask_ch: list[float] = []
        trmask_strips: list[float] = []

        trmask_string = adcdata_splitted[trmask_idx] if len(adcdata_splitted) > trmask_idx else ""
        trmask_vect_split = [token for token in trmask_string.split(trigger_mask_separator) if token != ""]

        # The trigger mask string is expected to contain channel numbers separated by underscores. 
        # We split the string by underscores and convert each token to an integer channel number. 
        # If the conversion fails or if the channel number is negative, we skip that token. 
        # For valid channel numbers, we add them to the list of trigger mask channels 
        # and also find their corresponding strip indices (if they exist in the stripIndices list) to add to the trigger mask strips. 
        # This allows us to keep track of which channels and strips were involved in the trigger for this board.
        for token in trmask_vect_split:
            trmask_channel = int(_safe_float(token, default=-1.0))
            if trmask_channel < 0:
                continue

            trmask_ch.append(float(trmask_channel))

            try:
                # The strip index is determined by finding the position of the channel number in the stripIndices list. 
                # If the channel number is not found in stripIndices, a ValueError is raised, and we catch it to ignore that channel for strip indexing. 
                # This mimics the C++ behavior where missing channels are silently ignored for strip indices.
                trmask_strip = stripIndices.index(trmask_channel)
                trmask_strips.append(float(trmask_strip))
            except ValueError:
                # Keep C++ behavior: missing channel is silently ignored for strip index.
                pass

        boards_trigger_masks_channels.append(trmask_ch)
        boards_trigger_masks_strips.append(trmask_strips)
        # How many strips were involved
        trigger_mask_sizes.append(len(trmask_strips))

        all_adc: list[float] = []
        for n_ch in range(n_channels):
            adc_idx = adc_offset + base_idx + n_ch
            adc_value = _safe_float(adcdata_splitted[adc_idx]) if len(adcdata_splitted) > adc_idx else 0.0
            all_adc.append(adc_value)
        # We need to reorder the ADC values according to the stripIndices, which represent the mapping from channel numbers to strip indices. 
        # For each strip index in stripIndices, we take the corresponding ADC value from all_adc. 
        # If a strip index is out of bounds for the all_adc list, we use a default value of 0.0. 
        # This ensures that the boards list contains ADC values ordered by strip index,
        sorted_all_adc = [all_adc[idx] if 0 <= idx < len(all_adc) else 0.0 for idx in stripIndices]
        boards.append(sorted_all_adc)

    return EventStruct(
        boards=boards,
        TrMask_channels=boards_trigger_masks_channels,
        TrMask_strips=boards_trigger_masks_strips,
        TrMask_size=trigger_mask_sizes,
        timeExp=all_time_exp,
        timeStamp=time_stamp,
    )


# Pythonic alias
read_event = ReadEvent
