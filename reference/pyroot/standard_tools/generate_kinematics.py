#!/usr/bin/env python3
"""
通用A(a,b)B 型反应计算
根据给定质心系角度计算运动学/或者随机生成
"""

from nuclear_data import get_particle_mass
import ROOT
import numpy as np
from array import array
from math import pi, sqrt, cos, sin, acos
import time
from datetime import datetime

def get_reaction_parameters(projectile_mass, projectile_energy,
                           target_type, target_energy=0,
                           ejectile1_mass=None, ejectile1_excitation=0,
                           ejectile2_mass=None, ejectile2_excitation=0,
                           theta_cm=None):
    """
    获取反应参数，可根据给定质心系角度计算运动学
    
    参数:
    ----------
    projectile_mass : float or str
        入射粒子质量(MeV)或符号
    projectile_energy : float
        入射粒子总动能(MeV)
    target_type : float or str
        靶核质量(MeV)或符号
    target_energy : float, default=0
        靶核动能(MeV)
    ejectile1_mass : float or str, optional
        出射粒子1质量，None表示与入射粒子相同
    ejectile1_excitation : float, default=0
        出射粒子1激发能(MeV)
    ejectile2_mass : float or str, optional
        出射粒子2质量，None表示与靶核相同
    ejectile2_excitation : float, default=0
        出射粒子2激发能(MeV)
    theta_cm : float, optional
        指定的出射粒子1质心系角度(度)，None表示随机抽样
    
    返回:
    ----------
    dict : 包含计算参数的字典
    'm_proj': m_proj,
    'E_proj_kin': projectile_energy,
    'm_target': m_target,
    'E_target_kin': target_energy,
    'm_ej1': m_ej1,
    'm_ej1_ex': m_ej1_ex,
    'excitation1': ejectile1_excitation,
    'm_ej2': m_ej2,
    'm_ej2_ex': m_ej2_ex,
    'excitation2': ejectile2_excitation,
    'Q_value': Q_value,
    'theta_cm': theta_cm,  # 指定的质心系角度
    'proj_label': str(projectile_mass),
    'target_label': str(target_type),
    'ej1_label': ej1_label,
    'ej2_label': ej2_label,
    """
    # 获取质量
    m_proj = get_particle_mass(projectile_mass)
    m_target = get_particle_mass(target_type)
    
    # 处理出射粒子
    if ejectile1_mass is None:
        m_ej1 = m_proj
        ej1_label = str(projectile_mass)
    else:
        m_ej1 = get_particle_mass(ejectile1_mass)
        ej1_label = str(ejectile1_mass)
    
    if ejectile2_mass is None:
        m_ej2 = m_target
        ej2_label = str(target_type)
    else:
        m_ej2 = get_particle_mass(ejectile2_mass)
        ej2_label = str(ejectile2_mass)
    
    # 激发态质量
    m_ej1_ex = m_ej1 + ejectile1_excitation
    m_ej2_ex = m_ej2 + ejectile2_excitation
    
    # Q值计算
    Q_value = (m_proj + m_target) - (m_ej1_ex + m_ej2_ex)
    
    return {
        'm_proj': m_proj,
        'E_proj_kin': projectile_energy,
        'm_target': m_target,
        'E_target_kin': target_energy,
        'm_ej1': m_ej1,
        'm_ej1_ex': m_ej1_ex,
        'excitation1': ejectile1_excitation,
        'm_ej2': m_ej2,
        'm_ej2_ex': m_ej2_ex,
        'excitation2': ejectile2_excitation,
        'Q_value': Q_value,
        'theta_cm': theta_cm,  # 指定的质心系角度
        'proj_label': str(projectile_mass),
        'target_label': str(target_type),
        'ej1_label': ej1_label,
        'ej2_label': ej2_label,
    }

def calculate_kinematics(params):
    """
    A(a,b)B性反应计算
    根据给定参数计算运动学（确定性的运动学计算）
    需要指定重核的质心系角度
    参数为get_reaction_parameters函数的返回值
    
    返回:
    ----------
    tuple : (theta_ej1_lab, E_ej1_lab, theta_ej2_lab, E_ej2_lab, theta_ej1_cm, E_ej1_cm)
          如果运动学不允许返回None
    """
    try:
        # 检查是否指定了质心系角度
        if params['theta_cm'] is None:
            raise ValueError("必须指定质心系角度(theta_cm)")
        
        theta_cm_rad = params['theta_cm'] * pi / 180.0  # 转换为弧度
        
        # 1. 计算初始总四矢量
        E_proj_total = params['E_proj_kin'] + params['m_proj']
        pz_proj = sqrt(E_proj_total**2 - params['m_proj']**2)
        beam = ROOT.TLorentzVector(0.0, 0.0, pz_proj, E_proj_total)
        
        # 2. 靶核四矢量
        if params['E_target_kin'] > 0:
            E_target_total = params['E_target_kin'] + params['m_target']
            pz_target = -sqrt(E_target_total**2 - params['m_target']**2)
            target = ROOT.TLorentzVector(0.0, 0.0, pz_target, E_target_total)
        else:
            target = ROOT.TLorentzVector(0.0, 0.0, 0.0, params['m_target'])
        
        # 3. 总四矢量和质心系
        total = beam + target
        sqrt_s = total.M()  # 质心系总能量
        
        # 4. 检查运动学可行性
        if sqrt_s < (params['m_ej1_ex'] + params['m_ej2_ex']):
            return None
        
        # 5. 计算质心系动量（使用相对论运动学公式）
        # 质心系总动量大小
        s = sqrt_s * sqrt_s  # Mandelstam s
        
        m3 = params['m_ej1_ex']  # 出射粒子1质量（激发态）
        m4 = params['m_ej2_ex']  # 出射粒子2质量（激发态）
        
        # 质心系动量公式: p_cm = sqrt([s - (m3+m4)^2][s - (m3-m4)^2]) / (2√s)
        p_cm = sqrt((s - (m3 + m4)**2) * (s - (m3 - m4)**2)) / (2.0 * sqrt_s)
        
        if p_cm <= 0:
            return None
        
        # 6. 计算质心系能量
        E3_cm = sqrt(p_cm**2 + m3**2)  # 出射粒子1质心系总能量
        E4_cm = sqrt(p_cm**2 + m4**2)  # 出射粒子2质心系总能量
        
        # 7. 创建质心系四矢量
        p3x_cm = p_cm * sin(theta_cm_rad)  # 假设在x-z平面
        p3z_cm = p_cm * cos(theta_cm_rad)
        
        ej1_cm = ROOT.TLorentzVector(p3x_cm, 0.0, p3z_cm, E3_cm)
        
        # 根据动量守恒，粒子2在质心系中方向相反
        ej2_cm = ROOT.TLorentzVector(-p3x_cm, 0.0, -p3z_cm, E4_cm)
        
        # 8. 变换到实验室系
        beta_cm = total.BoostVector()
        
        ej1_lab = ej1_cm.Clone()
        ej1_lab.Boost(beta_cm)
        
        ej2_lab = ej2_cm.Clone()
        ej2_lab.Boost(beta_cm)
        
        # 9. 计算物理量
        theta_ej1_lab = ej1_lab.Theta() * 180.0 / pi
        E_ej1_lab = ej1_lab.E() - m3  # 实验室系动能
        
        theta_ej2_lab = ej2_lab.Theta() * 180.0 / pi
        E_ej2_lab = ej2_lab.E() - m4  # 实验室系动能
        
        theta_ej1_cm = params['theta_cm']  # 给定的质心系角度
        E_ej1_cm = E3_cm - m3  # 质心系动能
        
        return (theta_ej1_lab, E_ej1_lab, 
                theta_ej2_lab, E_ej2_lab,
                theta_ej1_cm, E_ej1_cm)
        
    except Exception as e:
        print(f"运动学计算错误: {e}")
        return None
def simulate_single_event(params):
    """
    模拟单个反应事件 随机抽样
    参数为get_reaction_parameters函数的返回值
    tuple : (theta_ej1_lab, E_ej1_lab, theta_ej2_lab, E_ej2_lab, theta_ej1_cm, E_ej1_cm)
          如果运动学不允许返回None
    """
    # 如果指定了质心系角度，使用确定性的运动学计算
    if params.get('theta_cm') is not None:
        return calculate_kinematics(params)
    
    # 否则使用随机相空间抽样（保持向后兼容）
    try:
        E_proj_total = params['E_proj_kin'] + params['m_proj']
        pz_proj = ROOT.TMath.Sqrt(E_proj_total**2 - params['m_proj']**2)
        beam = ROOT.TLorentzVector(0.0, 0.0, pz_proj, E_proj_total)
        
        if params['E_target_kin'] > 0:
            E_target_total = params['E_target_kin'] + params['m_target']
            pz_target = -ROOT.TMath.Sqrt(E_target_total**2 - params['m_target']**2)
            target = ROOT.TLorentzVector(0.0, 0.0, pz_target, E_target_total)
        else:
            target = ROOT.TLorentzVector(0.0, 0.0, 0.0, params['m_target'])
        
        total = beam + target
        
        if total.M() < (params['m_ej1_ex'] + params['m_ej2_ex']):
            return None
        
        masses = array('d', [params['m_ej1_ex'], params['m_ej2_ex']])
        phase_space = ROOT.TGenPhaseSpace()
        
        if not phase_space.SetDecay(total, 2, masses):
            return None
        
        seed = int(time.time_ns() % 1000000)
        ROOT.gRandom.SetSeed(seed)
        
        weight = phase_space.Generate()
        if weight <= 0:
            return None
        
        ej1 = phase_space.GetDecay(0)
        ej2 = phase_space.GetDecay(1)
        
        theta_ej1_lab = ej1.Theta() * 180.0 / pi
        E_ej1_lab = ej1.E() - params['m_ej1_ex']
        
        theta_ej2_lab = ej2.Theta() * 180.0 / pi
        E_ej2_lab = ej2.E() - params['m_ej2_ex']
        
        beta_cm = total.BoostVector()
        ej1_cm = ej1.Clone()
        ej1_cm.Boost(-beta_cm)
        
        theta_ej1_cm = ej1_cm.Theta() * 180.0 / pi
        E_ej1_cm = ej1_cm.E() - params['m_ej1_ex']
        
        return (theta_ej1_lab, E_ej1_lab, 
                theta_ej2_lab, E_ej2_lab,
                theta_ej1_cm, E_ej1_cm)
        
    except Exception as e:
        print(f"模拟错误: {e}")
        return None


if __name__ == '__main__':
    # 示例用法
    parameters = get_reaction_parameters(
        projectile_mass='14O',
        projectile_energy=14*35,
        target_type='2H', # 靶核
        theta_cm=90
    )
    for i in parameters:
        print(i, parameters[i])

    kinematics = calculate_kinematics(parameters)
    print(kinematics[1]+kinematics[3])