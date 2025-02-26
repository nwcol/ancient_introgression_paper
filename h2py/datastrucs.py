
import numpy as np


class DpStats(list):

    x = 0


def Dp_names(num_pops):

    stat_names = []
    for i in range(num_pops):
        for j in range(i, num_pops):
            stat_names.append(f"Dplus_{i}_{j}")
    return stat_names


def H_names(num_pops):

    stat_names = []
    for i in range(num_pops):
        for j in range(i, num_pops):
            stat_names.append(f"H_{i}_{j}")
    return stat_names

