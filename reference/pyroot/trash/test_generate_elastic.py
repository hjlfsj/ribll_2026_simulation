#!/usr/bin/env python3
"""
优化版通用核反应模拟器：A(a,b)B 型反应
根据给定质心系角度计算运动学
"""

import ROOT
import numpy as np
import time
from datetime import datetime
from array import array
from math import pi, sqrt, cos, sin, acos

# ====================== 原子质量数据库 ======================
# 931.49410242 MeV/c^2
# Mass Excess经过校对，与nndc数据一致，第一项不一定准确
ATOMIC_MASS_DB = {
    '1H': (938.783, 7.288971),
    '2H': (1875.61294257, 13.135722),
    '3H': (2808.92113298, 14.949810),
    '4H': (3728.4, 24.6),
    
    '3He': (2808.391607, 14.931213),
    '4He': (3727.3794066, 2.424915),
    '5He': (4668.0, 11.231),
    '6He': (5605.577, 17.592),
    '8He': (7472.0, 31.61),
    
    '6Li': (5601.518, 14.086),
    '7Li': (6533.834, 14.908),
    '8Li': (7471.636, 20.946),
    '9Li': (8401.762, 24.954),
    '11Li': (10270.0, 40.728),
    
    '7Be': (6534.184, 15.769),
    '9Be': (8392.748, 11.348),
    '10Be': (9326.104, 12.607),
    '11Be': (10263.0, 20.17),
    
    '10B': (9324.436, 12.051),
    '11B': (10252.542, 8.668),
    '12B': (11178.5, 13.37),

    '10C': (0, 15.699),
    '11C': (0, 10.649),
    '12C': (11174.862, 0.0),
    '13C': (12109.482, 3.125),
    '14C': (13040.988, 3.020),
    '15C': (0, 9.873),
    
    '13N': (12110.6, 5.345),
    '14N': (13043.2, 2.863),
    '15N': (13971.8, 0.101),
    
    '14O': (13036.0, 8.008),
    '15O': (13971.2, 2.855),
    '16O': (14895.079, -4.737),
    '17O': (15828.0, -0.809),
    '18O': (16763.0, -0.782),
    
    
}

def get_particle_mass(particle):
    # if isinstance(particle, (int, float)):
    #     return float(particle)
    
    # particle = str(particle).strip()
    
    # if particle in ATOMIC_MASS_DB:
    #     return ATOMIC_MASS_DB[particle][0]
    
    try:
        import re
        match = re.match(r'(\d+)([A-Za-z]+)', particle)
        if match:
            isotope = match.group(1) + match.group(2)
            if isotope in ATOMIC_MASS_DB:
                return ATOMIC_MASS_DB[isotope][1] + int(match.group(1)) * 931.494
            else:
                A = int(match.group(1))
                u = 931.494
                return A * u
    except:
        pass
    
    raise ValueError(f"未知粒子符号: {particle}")

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
    根据给定参数计算运动学（确定性的运动学计算）
    
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
    模拟单个反应事件（根据给定角度或随机抽样）
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

# ====================== 主函数 ======================
def main():
    print("="*70)
    print("核反应模拟器 - TGraph运动学曲线")
    print("="*70)
    
    ROOT.gROOT.SetBatch(True)  # 批处理模式，不显示窗口
    
    excit_e = 0
    particle_A = '14O'
    particle_B = '2H'
    particle_C = '14O'
    particle_D = '2H'
    # 基本反应参数
    base_params = get_reaction_parameters(
        projectile_mass='14O',
        projectile_energy=14*35,
        target_type='2H',
        ejectile1_mass= particle_C,
        ejectile1_excitation=excit_e,
        ejectile2_mass= particle_D
    )
    
    print(f"反应参数:")
    print(f"  入射粒子: {base_params['proj_label']}, E = {base_params['E_proj_kin']} MeV")
    print(f"  靶核: {base_params['target_label']}")
    print(f"  出射粒子1: {base_params['ej1_label']}, Ex = {base_params['excitation1']} MeV")
    print(f"  出射粒子2: {base_params['ej2_label']}")
    print(f"  Q值: {base_params['Q_value']:.3f} MeV")
    print(base_params['m_ej1'])
    print(base_params['m_ej2'])
    # 按质心系角度步长生成结果
    print(f"\n按质心系角度生成运动学结果...")
    
    # 定义角度范围和步长
    theta_cm_min = 0.1
    theta_cm_max = 180.0
    theta_cm_step = 0.2
    
    n_angles = int((theta_cm_max - theta_cm_min) / theta_cm_step) 
    print(f"  角度范围: {theta_cm_min} 到 {theta_cm_max} 度")
    print(f"  步长: {theta_cm_step} 度")
    print(f"  总角度数: {n_angles}")
    
    start_time = time.time()
    
    # 存储所有角度的结果
    all_results = []
    
    # 循环所有质心系角度
    for i in range(n_angles):
        theta_cm = theta_cm_min + i * theta_cm_step
        
        # 每20度显示一次进度
        if i % 20 == 0:
            elapsed = time.time() - start_time
            if elapsed > 0:
                rate = i / elapsed
                print(f"  进度: {i}/{n_angles} ({i/n_angles*100:.1f}%), "
                      f"当前角度: {theta_cm:.1f}°, 速率: {rate:.0f} 角度/秒")
        
        # 为每个角度创建参数
        params = get_reaction_parameters(
        projectile_mass='14O',
        projectile_energy=14*35,
        target_type='2H',
        ejectile1_mass= particle_C,
        ejectile1_excitation=excit_e,
        ejectile2_mass= particle_D,
        theta_cm=theta_cm
    )
        
        # 计算运动学
        result = simulate_single_event(params)
        
        if result is not None:
            all_results.append(result)
    
    total_time = time.time() - start_time
    n_valid = len(all_results)
    
    print(f"\n计算完成!")
    print(f"总时间: {total_time:.2f} 秒")
    print(f"有效角度: {n_valid}/{n_angles} ({n_valid/n_angles*100:.1f}%)")
    print(f"平均速率: {n_angles/total_time:.1f} 角度/秒")
    
    # 准备数据用于TGraph
    n_points = n_valid
    
    # 提取数据
    theta_ej1_lab_list = [r[0] for r in all_results]
    E_ej1_lab_list = [r[1] for r in all_results]
    theta_ej2_lab_list = [r[2] for r in all_results]
    E_ej2_lab_list = [r[3] for r in all_results]

    theta_ej1_cm_list = [r[4] for r in all_results]
    E_ej1_cm_list = [r[5] for r in all_results]
    
    # 创建TGraph对象
    # 1. 出射粒子1实验室系能量-角度
    gr1 = ROOT.TGraph(n_points)
    gr1.SetTitle(f"{base_params['ej1_label']}* : E vs theta_lab  ;theta_lab [deg];E [MeV]")
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(0.5)
    gr1.SetMarkerColor(ROOT.kRed)
    gr1.SetLineColor(ROOT.kRed)
    gr1.SetLineWidth(2)
    
    # 2. 出射粒子2实验室系能量-角度
    gr2 = ROOT.TGraph(n_points)
    gr2.SetTitle(f"{base_params['ej2_label']} : E vs theta_lab  ;theta_lab [deg];E [MeV]")
    gr2.SetMarkerStyle(21)
    gr2.SetMarkerSize(0.5)
    gr2.SetMarkerColor(ROOT.kBlue)
    gr2.SetLineColor(ROOT.kBlue)
    gr2.SetLineWidth(2)
    
    # 3. 双粒子角度关联
    gr3 = ROOT.TGraph(n_points)
    gr3.SetTitle(f"Angular Correlation;theta_lab({base_params['ej1_label']}*) [deg];theta_lab({base_params['ej2_label']}) [deg]")
    gr3.SetMarkerStyle(22)
    gr3.SetMarkerSize(0.5)
    gr3.SetMarkerColor(ROOT.kGreen+2)
    gr3.SetLineColor(ROOT.kGreen+2)
    gr3.SetLineWidth(2)
    
    # 4. 出射粒子1质心系角度-实验室系能量
    gr4 = ROOT.TGraph(n_points)
    gr4.SetTitle(f"{base_params['ej1_label']}* CM Angle vs Lab Energy;theta_CM [deg];E_lab [MeV]")
    gr4.SetMarkerStyle(23)
    gr4.SetMarkerSize(0.5)
    gr4.SetMarkerColor(ROOT.kMagenta)
    gr4.SetLineColor(ROOT.kMagenta)
    gr4.SetLineWidth(2)
    
    # 5. 双粒子能量关联
    gr5 = ROOT.TGraph(n_points)
    gr5.SetTitle(f"Energy Correlation;E_lab({base_params['ej1_label']}*) [MeV];E_lab({base_params['ej2_label']}) [MeV]")
    gr5.SetMarkerStyle(24)
    gr5.SetMarkerSize(0.5)
    gr5.SetMarkerColor(ROOT.kOrange+1)
    gr5.SetLineColor(ROOT.kOrange+1)
    gr5.SetLineWidth(2)
    
    # 6. 实验室角度 vs 质心角度
    gr6 = ROOT.TGraph(n_points)
    gr6.SetTitle(f"Lab Angle vs CM Angle;θ_CM [deg];θ_lab [deg]")
    gr6.SetMarkerStyle(25)
    gr6.SetMarkerSize(0.5)
    gr6.SetMarkerColor(ROOT.kCyan+1)
    gr6.SetLineColor(ROOT.kCyan+1)
    gr6.SetLineWidth(2)
    
    # 填充TGraph数据
    for i, result in enumerate(all_results):
        (theta_ej1_lab, E_ej1_lab,
         theta_ej2_lab, E_ej2_lab,
         theta_ej1_cm, E_ej1_cm) = result
        
        gr1.SetPoint(i, theta_ej1_lab, E_ej1_lab)
        gr2.SetPoint(i, theta_ej2_lab, E_ej2_lab)
        gr3.SetPoint(i, theta_ej1_lab, theta_ej2_lab)
        gr4.SetPoint(i, theta_ej1_cm, E_ej1_lab)
        gr5.SetPoint(i, E_ej2_lab, E_ej1_lab)
        gr6.SetPoint(i, theta_ej1_cm, theta_ej1_lab)
    
    # 创建ROOT画布并绘制
    canvas = ROOT.TCanvas("canvas", "Kinematic Curves", 1600, 1200)
    canvas.Divide(3, 2)
    
    # 设置全局样式
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)
    
    # 绘制图1
    canvas.cd(1)
    gr1.Draw("AP")
    canvas.cd(1).Update()
    
    # 绘制图2
    canvas.cd(2)
    gr2.Draw("AP")
    canvas.cd(2).Update()
    
    # 绘制图3
    canvas.cd(3)
    gr3.Draw("AP")
    canvas.cd(3).Update()
    
    # 绘制图4
    canvas.cd(4)
    gr4.Draw("AP")
    canvas.cd(4).Update()
    
    # 绘制图5
    canvas.cd(5)
    gr5.Draw("AP")
    canvas.cd(5).Update()
    
    # 绘制图6
    canvas.cd(6)
    gr6.Draw("AP")
    canvas.cd(6).Update()
    
    canvas.Update()

    canvas.SaveAs("kinematic_curves.pdf")
    
    
    
    
    
    
    
 
    print("\n" + "="*70)
    print("运动学扫描统计:")
    print("="*70)
    
    if n_valid > 0:
        # 计算统计量
        theta_ej1_lab_arr = np.array(theta_ej1_lab_list)
        E_ej1_lab_arr = np.array(E_ej1_lab_list)
        theta_ej2_lab_arr = np.array(theta_ej2_lab_list)
        E_ej2_lab_arr = np.array(E_ej2_lab_list)
        
        print(f"\n{base_params['ej1_label']}* 实验室系:")
        print(f"  角度范围: {theta_ej1_lab_arr.min():.2f} 到 {theta_ej1_lab_arr.max():.2f} 度")
        print(f"  平均角度: {theta_ej1_lab_arr.mean():.1f} ± {theta_ej1_lab_arr.std():.1f} 度")
        print(f"  能量范围: {E_ej1_lab_arr.min():.2f} 到 {E_ej1_lab_arr.max():.2f} MeV")
        print(f"  平均能量: {E_ej1_lab_arr.mean():.1f} ± {E_ej1_lab_arr.std():.1f} MeV")
        
        print(f"\n{base_params['ej2_label']} 实验室系:")
        print(f"  角度范围: {theta_ej2_lab_arr.min():.1f} 到 {theta_ej2_lab_arr.max():.1f} 度")
        print(f"  平均角度: {theta_ej2_lab_arr.mean():.1f} ± {theta_ej2_lab_arr.std():.1f} 度")
        print(f"  能量范围: {E_ej2_lab_arr.min():.1f} 到 {E_ej2_lab_arr.max():.1f} MeV")
        print(f"  平均能量: {E_ej2_lab_arr.mean():.1f} ± {E_ej2_lab_arr.std():.1f} MeV")
        
        print(f"\n角度关联:")
        print(f"  当 θ_CM = 0° 时: θ_lab({base_params['ej1_label']}*) = {theta_ej1_lab_arr[0]:.1f}°, "
              f"θ_lab({base_params['ej2_label']}) = {theta_ej2_lab_arr[0]:.1f}°")
        print(f"  当 θ_CM = 90° 时: θ_lab({base_params['ej1_label']}*) = {theta_ej1_lab_arr[90]:.1f}°, "
              f"θ_lab({base_params['ej2_label']}) = {theta_ej2_lab_arr[90]:.1f}°")
        print(f"  当 θ_CM = 180° 时: θ_lab({base_params['ej1_label']}*) = {theta_ej1_lab_arr[-1]:.1f}°, "
              f"θ_lab({base_params['ej2_label']}) = {theta_ej2_lab_arr[-1]:.1f}°")
    
    print("\n" + "="*70)
    print("程序完成")

# ====================== 直接运行 ======================
if __name__ == "__main__":
    main()