from generate_kinematics import get_particle_mass
import ROOT
from math import sqrt, cos, sin, pi


# 缺失质量法A:  已知一个粒子的出射能量和角度，已知束流能量
def Q_missing_mass_mode1(A_name, a_name, b_name, B_name,
                   E_A_kin, E_b_kin, theta_b_deg):
    """
    使用四动量缺失质量法计算反应Q值
    
    参数:
    ----------
    A_name : str
        入射粒子名称 (如 '12C', '14O')
    a_name : str
        靶核名称 (如 '2H', '9Be')
    b_name : str
        测量粒子名称 (如 'p', '4He')
    B_name_gs : str
        剩余核的基态名称 (如 '13C', '15O')
    E_A_kin : float
        入射粒子A的动能 [MeV]
    E_b_kin : float
        出射粒子b的动能 [MeV]
    theta_b_deg : float
        出射粒子b与z轴的夹角 [度]
    
    返回:
    ----------
    tuple : (Q_measured, E_x_inferred, Q_gs_theoretical)
        Q_measured : 从运动学测量的Q值 [MeV]
        E_x_inferred : 推断的剩余核激发能 [MeV]（如果未知）
        Q_gs_theoretical : 基态反应的预期Q值 [MeV]
    """
    # 1. 获取粒子基态质量
    m_A_gs = get_particle_mass(A_name)
    m_a_gs = get_particle_mass(a_name)
    m_b_gs = get_particle_mass(b_name)
    m_B_gs = get_particle_mass(B_name)
    
    # 2. 计算基态反应的预期Q值
    Q_gs_theoretical = (m_A_gs + m_a_gs) - (m_b_gs + m_B_gs)
    
    # 3. 创建四矢量
    # 入射粒子A
    E_A_total = E_A_kin + m_A_gs
    pz_A = sqrt(E_A_total**2 - m_A_gs**2)
    P_A = ROOT.TLorentzVector(0.0, 0.0, pz_A, E_A_total)
    
    # 靶核a（静止）
    P_a = ROOT.TLorentzVector(0.0, 0.0, 0.0, m_a_gs)
    
    # 出射粒子b
    theta_b_rad = theta_b_deg * pi / 180.0
    E_b_total = E_b_kin + m_b_gs
    p_b = sqrt(E_b_total**2 - m_b_gs**2)
    pz_b = p_b * cos(theta_b_rad)
    px_b = p_b * sin(theta_b_rad)
    P_b = ROOT.TLorentzVector(px_b, 0.0, pz_b, E_b_total)
    
    # 4. 计算剩余核B的四矢量
    P_B = P_A + P_a - P_b
    
    # 5. 测量剩余核的质量
    M_B_measured = P_B.M()
    
    # 6. 计算测量的Q值
    Q_measured = (m_A_gs + m_a_gs) - (m_b_gs + M_B_measured)
    
    # 7. 推断激发能
    # 如果剩余核处于激发态，其质量 = m_B_gs + E_x
    # 所以：M_B_measured = m_B_gs + E_x
    # 因此：E_x = M_B_measured - m_B_gs
    E_x_inferred = M_B_measured - m_B_gs
    
    
    
    # 9. 输出详细信息
    
    
    return Q_measured, E_x_inferred, Q_gs_theoretical





# test
if __name__ == "__main__":
    from generate_kinematics import get_reaction_parameters
    from generate_kinematics import simulate_single_event
    # 参数设置
    parameters = get_reaction_parameters(
        projectile_mass='14O',
        projectile_energy=14*35,
        target_type='2H', # 靶核
        ejectile1_mass='10C', 
        ejectile1_excitation=0,
        ejectile2_mass= '6Li',  
        ejectile2_excitation=0,
        theta_cm=9
    
        
    )
    for i in parameters:
        print(i,":  ", parameters[i])
    kinematics = simulate_single_event(parameters)
    # 注意返回，和参数设置保持一致
    # ----------
    # tuple : (theta_ej1_lab, E_ej1_lab, theta_ej2_lab, E_ej2_lab, theta_ej1_cm, E_ej1_cm)
    #       如果运动学不允许返回None
    
    print ('='*70)
    E_total = 14*35
    theta_b = kinematics[2]
    E_b = kinematics[3]

    Q_measured, E_x_inferred, Q_gs_theoretical = Q_missing_mass_mode1('14O', '2H', '6Li', '10C',
                        E_A_kin=E_total, E_b_kin=E_b, theta_b_deg=theta_b)
    print(f"Q_measured: {Q_measured:.3f} MeV")
    print(f"E_x_inferred: {E_x_inferred:.3f} MeV")
    print(f"Q_gs_theoretical: {Q_gs_theoretical:.3f} MeV")