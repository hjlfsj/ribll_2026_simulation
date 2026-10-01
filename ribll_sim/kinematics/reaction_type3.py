"""
反应类型3：A + B -> C* + D,  C* -> E + F
先进行类型2反应生成激发态C*和D，然后C*衰变成E和F。
E默认基态，F可设置激发态。

关键：E,F的方向计算在C*静止系中进行，以C*方向为z轴，
然后转换到实验室系（target的z轴）。
"""

import ROOT
from array import array
from math import pi, sqrt
import time

from ..data.nuclear_data import get_particle_mass
from .reaction_type2 import simulate_type2_reaction


def simulate_type3_sequential(particle_A, E_A_kin, particle_B,
                              particle_C, particle_D,
                              particle_E, particle_F,
                              excitation_C=0.0,
                              excitation_F=0.0,
                              cms_theta=180.0):
    """
    A + B -> C* + D,  C* -> E + F 级联反应模拟

    参数:
    ----------
    particle_A : str
        入射粒子符号
    E_A_kin : float
        入射粒子动能 [MeV]
    particle_B : str
        靶核符号
    particle_C : str
        中间产物C符号 (激发态)
    particle_D : str
        出射粒子D符号 (默认基态)
    particle_E : str
        C衰变产物E符号 (默认基态)
    particle_F : str
        C衰变产物F符号 (可设置激发态)
    excitation_C : float
        C的激发能 [MeV]
    excitation_F : float
        F的激发能 [MeV], 默认0
    cms_theta : float
        第一步产物C*在质心系中的最大极角 [度], (0, cms_theta) 均匀抽样

    返回:
    ----------
    dict : {
        'valid': bool,
        'step1': 类型2反应结果,
        'E': 产物E信息,
        'F': 产物F信息,
        'D': 产物D信息,
    }
    """
    step1 = simulate_type2_reaction(
        particle_A, E_A_kin, particle_B,
        particle_C, particle_D,
        excitation_C=excitation_C,
        cms_theta=cms_theta
    )

    if not step1['valid']:
        return {
            'valid': False,
            'reason': '第一步反应失败: ' + step1.get('reason', '未知'),
            'step1': step1
        }

    p4_Cstar = step1['C']['p4']
    p4_D = step1['D']['p4']

    m_E_gs = get_particle_mass(particle_E)
    m_F_gs = get_particle_mass(particle_F)

    m_Cstar = p4_Cstar.M()
    m_E = m_E_gs
    m_F = m_F_gs + excitation_F

    Q_breakup = m_Cstar - (m_E + m_F)

    if Q_breakup < 0:
        return {
            'valid': False,
            'reason': f'C*衰变Q值负值: {Q_breakup:.3f} MeV',
            'step1': step1
        }

    masses = array('d', [m_E, m_F])
    phase_space = ROOT.TGenPhaseSpace()

    if not phase_space.SetDecay(p4_Cstar, 2, masses):
        return {
            'valid': False,
            'reason': 'C*衰变相空间设置失败',
            'step1': step1
        }

    seed = int(time.time_ns() % 1000000)
    ROOT.gRandom.SetSeed(seed)
    weight = phase_space.Generate()

    if weight <= 0:
        return {
            'valid': False,
            'reason': 'C*衰变相空间生成失败',
            'step1': step1
        }

    p4_E_lab = phase_space.GetDecay(0)
    p4_F_lab = phase_space.GetDecay(1)

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
        'step1': step1,
        'cms_theta_C': step1.get('cms_theta_C', None),
        'Q_breakup': Q_breakup,
        'E': make_particle_info(particle_E, p4_E_lab, m_E_gs, 0.0),
        'F': make_particle_info(particle_F, p4_F_lab, m_F_gs, excitation_F),
        'D': {
            'name': particle_D,
            'mass_gs': step1['D']['mass_gs'],
            'excitation': 0.0,
            'mass': step1['D']['mass'],
            'p4': ROOT.TLorentzVector(p4_D),
            'Ek': step1['D']['Ek'],
            'theta': step1['D']['theta'],
            'phi': step1['D']['phi'],
            'p': step1['D']['p'],
        },
    }