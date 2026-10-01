"""
反应类型1：A* -> B + C
母核A处于激发态，衰变成B和C两个粒子。
A具有初始动能（沿z轴），B和C均为基态。
激发能加在母核A上，且必须大于衰变基态Q值的绝对值（-Q_gs）。
"""

import ROOT
from array import array
from math import pi, sqrt
import time
import random

from ..data.nuclear_data import get_particle_mass


def simulate_type1_decay(particle_A, E_A_kin, particle_B, particle_C,
                         excitation_C=0.0):
    """
    A* -> B + C 衰变模拟

    excitation_C 在此表示母核A的激发能。
    B和C始终为基态，A的质量 = 基态质量 + excitation_C。
    要求 excitation_C > -Q_gs（基态Q值的绝对值）。

    参数:
    ----------
    particle_A : str
        母核符号 (如 '14O')
    E_A_kin : float
        母核动能 [MeV]
    particle_B : str
        衰变产物B符号 (始终基态)
    particle_C : str
        衰变产物C符号 (始终基态)
    excitation_C : float
        母核A的激发能 [MeV], 默认0

    返回:
    ----------
    dict : {
        'valid': bool,
        'A': 母核信息,
        'B': 产物B信息,
        'C': 产物C信息,
        'Q_value': float
    }
    """
    m_A_gs = get_particle_mass(particle_A)
    m_B_gs = get_particle_mass(particle_B)
    m_C_gs = get_particle_mass(particle_C)

    m_B = m_B_gs
    m_C = m_C_gs

    m_A = m_A_gs + excitation_C

    Q_gs = m_A_gs - m_B_gs - m_C_gs
    Q_value = m_A - m_B - m_C

    if Q_gs < 0 and excitation_C <= -Q_gs:
        return {
            'valid': False,
            'reason': f'激发能不足以克服负Q值: 需要 >{-Q_gs:.3f} MeV, 实际 {excitation_C:.3f} MeV',
            'Q_value': Q_value
        }

    if Q_value <= 0:
        return {
            'valid': False,
            'reason': f'Q值非正: {Q_value:.3f} MeV (激发能={excitation_C:.3f} MeV)',
            'Q_value': Q_value
        }

    E_A_total = E_A_kin + m_A
    p_A = sqrt(E_A_total**2 - m_A**2)
    P_A = ROOT.TLorentzVector(0.0, 0.0, p_A, E_A_total)

    masses = array('d', [m_B, m_C])
    phase_space = ROOT.TGenPhaseSpace()

    if not phase_space.SetDecay(P_A, 2, masses):
        return {
            'valid': False,
            'reason': '相空间设置失败',
            'Q_value': Q_value
        }

    seed = int(time.time_ns() % 1000000)
    ROOT.gRandom.SetSeed(seed)
    weight = phase_space.Generate()

    if weight <= 0:
        return {
            'valid': False,
            'reason': '相空间生成失败',
            'Q_value': Q_value
        }

    p4_B = phase_space.GetDecay(0)
    p4_C = phase_space.GetDecay(1)

    def make_particle_info(name, p4, mass_gs, excitation):
        return {
            'name': name,
            'mass_gs': mass_gs,
            'excitation': excitation,
            'mass': mass_gs + excitation,
            'p4': ROOT.TLorentzVector(p4),
            'Ek': p4.E() - (mass_gs + excitation),
            'theta': p4.Theta() * 180.0 / pi,
            'phi': p4.Phi() * 180.0 / pi,
            'p': p4.P(),
        }

    return {
        'valid': True,
        'Q_value': Q_value,
        'A': {
            'name': particle_A,
            'mass_gs': m_A_gs,
            'excitation': excitation_C,
            'mass': m_A,
            'Ek': E_A_kin,
            'p4': ROOT.TLorentzVector(P_A),
        },
        'B': make_particle_info(particle_B, p4_B, m_B_gs, 0.0),
        'C': make_particle_info(particle_C, p4_C, m_C_gs, 0.0),
    }