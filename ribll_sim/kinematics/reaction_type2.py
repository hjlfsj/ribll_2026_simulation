"""
反应类型2：A + B -> C + D
A以一定动能入射，与静止的B碰撞生成C和D。
A, B默认基态。C可设置激发态，D默认基态。
当C=A, D=B时为弹散/非弹散；否则为转移反应。

质心系采样：以产物C的质心系角度在(0, cms_theta)中均匀抽样。
"""

import ROOT
from math import pi, sqrt, sin, cos
import random

from ..data.nuclear_data import get_particle_mass


def simulate_type2_reaction(particle_A, E_A_kin, particle_B,
                            particle_C, particle_D,
                            excitation_C=0.0,
                            cms_theta=180.0):
    """
    A + B -> C + D 反应模拟

    参数:
    ----------
    particle_A : str
        入射粒子符号 (如 '14O')
    E_A_kin : float
        入射粒子动能 [MeV]
    particle_B : str
        靶核符号 (如 '2H')
    particle_C : str
        出射粒子C符号 (可设置激发态)
    particle_D : str
        出射粒子D符号 (默认基态)
    excitation_C : float
        C的激发能 [MeV], 默认0
    cms_theta : float
        产物C在质心系中的最大极角 [度], (0, cms_theta) 均匀抽样, 默认180

    返回:
    ----------
    dict : {
        'valid': bool,
        'C': 产物C信息,
        'D': 产物D信息,
        'Q_value': float,
        'reaction_type': 'elastic' | 'inelastic' | 'transfer'
    }
    """
    m_A = get_particle_mass(particle_A)
    m_B = get_particle_mass(particle_B)
    m_C_gs = get_particle_mass(particle_C)
    m_D_gs = get_particle_mass(particle_D)

    m_C = m_C_gs + excitation_C
    m_D = m_D_gs

    Q_value = (m_A + m_B) - (m_C + m_D)

    if (particle_A == particle_C) and (particle_B == particle_D):
        if excitation_C == 0:
            reaction_type = 'elastic'
        else:
            reaction_type = 'inelastic'
    else:
        reaction_type = 'transfer'

    if Q_value < 0:
        required_E = -Q_value * (m_A + m_B) / m_B
        if E_A_kin < required_E:
            return {
                'valid': False,
                'reason': f'动能不足以克服负Q值: 需要>{required_E:.1f} MeV, 实际{E_A_kin:.1f} MeV',
                'Q_value': Q_value,
                'reaction_type': reaction_type
            }

    E_A_total = E_A_kin + m_A
    p_A = sqrt(E_A_total**2 - m_A**2)
    P_A = ROOT.TLorentzVector(0.0, 0.0, p_A, E_A_total)
    P_B = ROOT.TLorentzVector(0.0, 0.0, 0.0, m_B)

    P_total = P_A + P_B
    sqrt_s = P_total.M()

    if sqrt_s < (m_C + m_D):
        return {
            'valid': False,
            'reason': '质心系能量不足以生成产物',
            'Q_value': Q_value,
            'reaction_type': reaction_type
        }

    beta_cm = P_total.BoostVector()

    s = sqrt_s * sqrt_s
    p_cm = sqrt(max(0, (s - (m_C + m_D)**2) * (s - (m_C - m_D)**2))) / (2.0 * sqrt_s)

    if p_cm <= 0:
        return {
            'valid': False,
            'reason': '质心系动量计算失败',
            'Q_value': Q_value,
            'reaction_type': reaction_type
        }

    E_C_cm = sqrt(p_cm**2 + m_C**2)
    E_D_cm = sqrt(p_cm**2 + m_D**2)

    theta_cm_rad = random.uniform(0, cms_theta) * pi / 180.0
    phi_cm = random.uniform(0, 2.0 * pi)

    px_cm = p_cm * sin(theta_cm_rad) * cos(phi_cm)
    py_cm = p_cm * sin(theta_cm_rad) * sin(phi_cm)
    pz_cm = p_cm * cos(theta_cm_rad)

    p4_C_cm = ROOT.TLorentzVector(px_cm, py_cm, pz_cm, E_C_cm)
    p4_D_cm = ROOT.TLorentzVector(-px_cm, -py_cm, -pz_cm, E_D_cm)

    p4_C_lab = p4_C_cm.Clone()
    p4_C_lab.Boost(beta_cm)
    p4_D_lab = p4_D_cm.Clone()
    p4_D_lab.Boost(beta_cm)

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
        'reaction_type': reaction_type,
        'cms_theta_C': theta_cm_rad * 180.0 / pi,
        'C': make_particle_info(particle_C, p4_C_lab, m_C_gs, excitation_C),
        'D': make_particle_info(particle_D, p4_D_lab, m_D_gs, 0.0),
    }