"""
激发能重建模块
从测量的粒子能量和角度重建反应激发能谱
"""

from math import sqrt, pi, sin, cos
import ROOT

from ..data.nuclear_data import get_particle_mass


def reconstruct_excitation_type1(result, particle_A, particle_B, particle_C):
    """
    类型1: 重建母核A的激发能
    A -> B + C, 测量B和C的能量和角度

    E_x(A) = M_reconstructed(A) - m_A_gs
    """
    m_A_gs = get_particle_mass(particle_A)

    p4_B = result['B']['p4']
    p4_C = result['C']['p4']
    p4_total = p4_B + p4_C

    M_reco = p4_total.M()
    E_x = M_reco - m_A_gs
    return max(0.0, E_x)


def reconstruct_excitation_type2(result, particle_A, E_beam, particle_B,
                                 particle_C, particle_D):
    """
    类型2: 重建C的激发能
    A + B -> C + D, 测量C和D

    使用缺失质量法: P_C = P_A + P_B - P_D
    E_x(C) = M_C_reconstructed - m_C_gs
    """
    m_A = get_particle_mass(particle_A)
    m_B = get_particle_mass(particle_B)
    m_C_gs = get_particle_mass(particle_C)

    E_A_total = E_beam + m_A
    p_A = sqrt(E_A_total**2 - m_A**2)
    P_A = ROOT.TLorentzVector(0.0, 0.0, p_A, E_A_total)
    P_B = ROOT.TLorentzVector(0.0, 0.0, 0.0, m_B)

    p4_D = result['D']['p4']
    P_C = P_A + P_B - p4_D

    M_C_reco = P_C.M()
    E_x = M_C_reco - m_C_gs
    return max(0.0, E_x)


def reconstruct_excitation_type3(result, E_beam, particle_A, particle_B,
                                 particle_C, particle_D, particle_F):
    """
    类型3: 重建C的激发能
    A+B -> C*+D -> (E+F)+D

    使用不变质量法: P_Cstar = P_E + P_F
    E_x(C) = M_E+F - m_C_gs
    """
    m_C_gs = get_particle_mass(particle_C)

    p4_E = result['E']['p4']
    p4_F = result['F']['p4']
    P_Cstar = p4_E + p4_F

    M_Cstar = P_Cstar.M()
    E_x = M_Cstar - m_C_gs
    return max(0.0, E_x)


def reconstruct_excitation(kin_result, params):
    """根据反应类型自动选择重建方法"""
    rt = params['reaction_type']
    try:
        if rt == 1:
            return reconstruct_excitation_type1(
                kin_result,
                params['particle_A'], params['particle_B'], params['particle_C']
            )
        elif rt == 2:
            return reconstruct_excitation_type2(
                kin_result,
                params['particle_A'], params['E_beam'],
                params['particle_B'], params['particle_C'], params['particle_D']
            )
        elif rt == 3:
            return reconstruct_excitation_type3(
                kin_result,
                params['E_beam'],
                params['particle_A'], params['particle_B'],
                params['particle_C'], params['particle_D'], params['particle_F']
            )
    except Exception:
        return 0.0
    return 0.0


def _build_lorentz_from_ek_theta_phi(Ek, theta_rad, phi_rad, particle_name):
    """从动能和角度构建TLorentzVector"""
    m = get_particle_mass(particle_name)
    E_tot = Ek + m
    p_val = sqrt(max(0.0, E_tot**2 - m**2))
    return ROOT.TLorentzVector(
        p_val * sin(theta_rad) * cos(phi_rad),
        p_val * sin(theta_rad) * sin(phi_rad),
        p_val * cos(theta_rad),
        E_tot,
    )


def reconstruct_excitation_experimental(detected_data, kin_result, params):
    """使用探测能量和角度重建激发能（含探测器分辨模糊效应）"""
    rt = params['reaction_type']
    m_A = get_particle_mass(params['particle_A'])
    m_B = get_particle_mass(params['particle_B'])
    E_beam = params['E_beam']

    E_A_total = E_beam + m_A
    p_A = sqrt(E_A_total**2 - m_A**2)
    P_A = ROOT.TLorentzVector(0.0, 0.0, p_A, E_A_total)
    P_B = ROOT.TLorentzVector(0.0, 0.0, 0.0, m_B)

    try:
        if rt == 1:
            m_A_gs = get_particle_mass(params['particle_A'])
            Ek_B = detected_data.get('det_p0_Eexp', detected_data.get('det_p0_Ek', 0))
            theta_B = detected_data.get('det_p0_theta', 0) * pi / 180.0
            phi_B = kin_result['B']['phi'] * pi / 180.0
            Ek_C = detected_data.get('det_p1_Eexp', detected_data.get('det_p1_Ek', 0))
            theta_C = detected_data.get('det_p1_theta', 0) * pi / 180.0
            phi_C = kin_result['C']['phi'] * pi / 180.0

            P_B_det = _build_lorentz_from_ek_theta_phi(
                Ek_B, theta_B, phi_B, kin_result['B']['name'])
            P_C_det = _build_lorentz_from_ek_theta_phi(
                Ek_C, theta_C, phi_C, kin_result['C']['name'])
            P_A_reco = P_B_det + P_C_det
            M_A_reco = P_A_reco.M()
            return max(0.0, M_A_reco - m_A_gs)

        elif rt == 2:
            m_C_gs = get_particle_mass(params['particle_C'])
            Ek_D = detected_data.get('det_p1_Eexp', detected_data.get('det_p1_Ek', 0))
            theta_D = detected_data.get('det_p1_theta', 0) * pi / 180.0
            phi_D = kin_result['D']['phi'] * pi / 180.0

            P_D_det = _build_lorentz_from_ek_theta_phi(
                Ek_D, theta_D, phi_D, params['particle_D'])
            P_C_reco = P_A + P_B - P_D_det
            M_C_reco = P_C_reco.M()
            return max(0.0, M_C_reco - m_C_gs)

        elif rt == 3:
            m_C_gs = get_particle_mass(params['particle_C'])
            Ek_E = detected_data.get('det_p0_Eexp', detected_data.get('det_p0_Ek', 0))
            theta_E = detected_data.get('det_p0_theta', 0) * pi / 180.0
            phi_E = kin_result['E']['phi'] * pi / 180.0
            Ek_F = detected_data.get('det_p1_Eexp', detected_data.get('det_p1_Ek', 0))
            theta_F = detected_data.get('det_p1_theta', 0) * pi / 180.0
            phi_F = kin_result['F']['phi'] * pi / 180.0

            P_E_det = _build_lorentz_from_ek_theta_phi(
                Ek_E, theta_E, phi_E, params['particle_E'])
            P_F_det = _build_lorentz_from_ek_theta_phi(
                Ek_F, theta_F, phi_F, params['particle_F'])
            P_Cstar = P_E_det + P_F_det
            M_Cstar = P_Cstar.M()
            return max(0.0, M_Cstar - m_C_gs)
    except Exception as e:
        import sys
        print("[ERROR] reconstruct_excitation_experimental:", e, file=sys.stderr)
        return 0.0
    return 0.0