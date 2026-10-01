#!/usr/bin/env python3
"""
多体衰变运动学模拟：激发态原子核衰变成多个粒子
"""
from nuclear_data import get_particle_mass
import ROOT
import numpy as np
from array import array
from math import pi, sqrt
import random

def breakup_kinematics(mass_A, Ek_A, theta_A, E_exciting_A, list_particle, list_ex, phi_A=0.0):
    """
    模拟激发态原子核的多体衰变
    
    参数:
    ----------
    mass_A : str or float
        母核种类或质量 [MeV]
    Ek_A : float
        母核动能 [MeV]
    theta_A : float
        母核运动方向与z轴夹角 [度]，范围: 0-180°
    E_exciting_A : float
        母核激发能 [MeV]
    list_particle : list of str
        衰变粒子种类列表，长度2-5
    list_ex : list of float
        对应衰变粒子的激发能列表 [MeV]
    phi_A : float, optional
        母核在x-y平面的方位角 [度]，默认0°，范围: 0-360°
        物理意义：从x轴正方向逆时针旋转的角度
    
    返回:
    ----------
    dict : 包含衰变粒子信息的字典，如果衰变不可能返回None
    """
    
    # 1. 参数验证
    n_particles = len(list_particle)
    if n_particles < 2 or n_particles > 5:
        raise ValueError(f"粒子数量必须为2-5，当前为 {n_particles}")
    
    if len(list_ex) != n_particles:
        raise ValueError(f"激发能列表长度必须与粒子列表相同")
    
    # 验证角度范围
    if not (0 <= theta_A <= 180):
        raise ValueError(f"极角theta_A必须在0-180度范围内，当前为 {theta_A}")
    
    # phi_A在0-360度范围内规范化
    phi_A = phi_A % 360.0  # 确保在0-360度范围内
    
    # 2. 获取母核质量
    if isinstance(mass_A, str):
        m_A_gs = get_particle_mass(mass_A)  # 基态质量
        A_name = mass_A
    else:
        m_A_gs = float(mass_A)
        A_name = f"母核({m_A_gs:.1f}MeV)"
    
    # 激发态质量
    m_A_ex = m_A_gs + E_exciting_A
    
    # 3. 获取衰变粒子质量
    m_daughters_gs = []
    m_daughters_ex = []
    
    for i, (particle, ex) in enumerate(zip(list_particle, list_ex)):
        m_gs = get_particle_mass(particle)
        m_ex = m_gs + ex
        m_daughters_gs.append(m_gs)
        m_daughters_ex.append(m_ex)
    
    # 4. 计算衰变Q值
    total_mass_daughters = sum(m_daughters_ex)
    Q_value = m_A_ex - total_mass_daughters
    
    if Q_value < 0:
        print(f"衰变不可能: Q值 = {Q_value:.3f} MeV < 0")
        print(f"母核质量: {m_A_ex:.3f} MeV")
        print(f"子核总质量: {total_mass_daughters:.3f} MeV")
        return {
            'valid': False,
            'reason': f'Q值负值: {Q_value:.3f} MeV',
            'Q_value': Q_value,
            'mother_mass': m_A_ex,
            'daughters_total_mass': total_mass_daughters
        }
    
    # 5. 创建母核在实验室系的四矢量（现在包含phi角）
    E_A_total = Ek_A + m_A_ex
    p_A = sqrt(E_A_total**2 - m_A_ex**2)
    
    # 将动量分解到三维空间
    theta_A_rad = theta_A * pi / 180.0
    phi_A_rad = phi_A * pi / 180.0
    
    # 球坐标到笛卡尔坐标转换：
    # x = p * sinθ * cosφ
    # y = p * sinθ * sinφ  
    # z = p * cosθ
    px_A = p_A * np.sin(theta_A_rad) * np.cos(phi_A_rad)
    py_A = p_A * np.sin(theta_A_rad) * np.sin(phi_A_rad)
    pz_A = p_A * np.cos(theta_A_rad)
    
    P_A = ROOT.TLorentzVector(px_A, py_A, pz_A, E_A_total)
    
    # 6. 使用TGenPhaseSpace生成多体衰变
    masses_array = array('d', m_daughters_ex)
    phase_space = ROOT.TGenPhaseSpace()
    
    if not phase_space.SetDecay(P_A, n_particles, masses_array):
        return {
            'valid': False,
            'reason': '相空间生成器设置失败',
            'Q_value': Q_value
        }
    
    # 7. 生成随机衰变事件
    seed = int(random.random() * 1000000)
    ROOT.gRandom.SetSeed(seed)
    
    weight = phase_space.Generate()
    if weight <= 0:
        return {
            'valid': False,
            'reason': '相空间生成失败',
            'Q_value': Q_value
        }
    
    # 8. 提取衰变粒子信息
    particles_info = []
    
    for i in range(n_particles):
        particle_vec = phase_space.GetDecay(i)
        
        # 计算实验室系角度
        theta_lab = particle_vec.Theta() * 180.0 / pi
        phi_lab = particle_vec.Phi() * 180.0 / pi
        
        # 确保phi_lab在0-360度范围内
        if phi_lab < 0:
            phi_lab += 360.0
        
        # 动能
        E_total = particle_vec.E()
        Ek_lab = E_total - m_daughters_ex[i]
        
        # 动量分量
        px = particle_vec.Px()
        py = particle_vec.Py()
        pz = particle_vec.Pz()
        
        particle_info = {
            'name': list_particle[i],
            'mass': m_daughters_gs[i],
            'excitation': list_ex[i],
            'mass_excited': m_daughters_ex[i],
            'theta_lab': theta_lab,
            'phi_lab': phi_lab,
            'Ek_lab': Ek_lab,
            'E_total': E_total,
            'px': px,
            'py': py,
            'pz': pz,
            'p': sqrt(px**2 + py**2 + pz**2),
            'index': i
        }
        
        particles_info.append(particle_info)
    
    # 9. 验证守恒
    total_p_initial = P_A
    total_p_final = ROOT.TLorentzVector(0, 0, 0, 0)
    
    for info in particles_info:
        p_vec = ROOT.TLorentzVector(info['px'], info['py'], info['pz'], info['E_total'])
        total_p_final += p_vec
    
    momentum_diff = (total_p_initial - total_p_final).P()
    energy_diff = abs(total_p_initial.E() - total_p_final.E())
    
    # 10. 返回结果（现在包含phi_A）
    result = {
        'valid': True,
        'mother_name': A_name,
        'mother_mass_gs': m_A_gs,
        'mother_mass_ex': m_A_ex,
        'mother_Ek': Ek_A,
        'mother_theta': theta_A,
        'mother_phi': phi_A,  # 新增
        'mother_excitation': E_exciting_A,
        'total_energy': E_A_total,
        'Q_value': Q_value,
        'n_particles': n_particles,
        'particles': particles_info,
        'weight': weight,
        'conservation_check': {
            'momentum_diff': momentum_diff,
            'energy_diff': energy_diff,
            'is_conserved': momentum_diff < 1e-6 and energy_diff < 1e-6
        }
    }
    
    return result

def print_decay_result(result):
    """格式化打印衰变结果"""
    if not result['valid']:
        print(f"衰变无效: {result.get('reason', '未知原因')}")
        return
    
    print("="*80)
    print("多体衰变模拟结果")
    print("="*80)
    
    print(f"母核: {result['mother_name']}")
    print(f"  基态质量: {result['mother_mass_gs']:.3f} MeV")
    print(f"  激发能: {result['mother_excitation']:.3f} MeV")
    print(f"  激发态质量: {result['mother_mass_ex']:.3f} MeV")
    print(f"  动能: {result['mother_Ek']:.3f} MeV, 角度: {result['mother_theta']:.1f}°")
    print(f"  总能量: {result['total_energy']:.3f} MeV")
    
    print(f"\n衰变Q值: {result['Q_value']:.3f} MeV")
    print(f"衰变粒子数: {result['n_particles']}")
    
    print(f"\n衰变粒子信息:")
    for i, particle in enumerate(result['particles']):
        print(f"  粒子{i+1}: {particle['name']}")
        print(f"    质量: {particle['mass']:.3f} MeV, 激发能: {particle['excitation']:.3f} MeV")
        print(f"    总质量: {particle['mass_excited']:.3f} MeV")
        print(f"    实验室系: Ek = {particle['Ek_lab']:.3f} MeV")
        print(f"    角度: θ = {particle['theta_lab']:.2f}°, φ = {particle['phi_lab']:.2f}°")
        print(f"    动量: p = {particle['p']:.3f} MeV/c")
        print(f"    动量分量: px = {particle['px']:.3f}, py = {particle['py']:.3f}, pz = {particle['pz']:.3f}")
        print()
    
    print(f"守恒检查:")
    check = result['conservation_check']
    print(f"  动量差异: {check['momentum_diff']:.6e} MeV/c")
    print(f"  能量差异: {check['energy_diff']:.6e} MeV")
    print(f"  守恒状态: {'OK' if check['is_conserved'] else '警告'}")

def generate_multiple_events(mass_A, Ek_A, theta_A, E_exciting_A, 
                           list_particle, list_ex,phi_A=0, n_events=100):
    """
    生成多个衰变事件
    
    返回:
    ----------
    list : 包含多个衰变事件的列表
    """
    events = []
    
    for event_id in range(n_events):
        result = breakup_kinematics(mass_A, Ek_A, theta_A, E_exciting_A, 
                                   list_particle, list_ex,phi_A)
        
        if result['valid']:
            result['event_id'] = event_id
            events.append(result)
    
    return events

def analyze_decay_distribution(events, particle_index=0):
    """
    分析衰变粒子的角度和能量分布
    """
    if not events:
        return None
    
    thetas = []
    energies = []
    
    new_thetas = []

    for event in events:
        particle = event['particles'][particle_index]
        thetas.append(particle['theta_lab'])
        energies.append(particle['Ek_lab'])
        new_thetas.append(transform_to_mother_frame(event)['particles'][particle_index]['theta_mother'])
    
    thetas = np.array(thetas)
    energies = np.array(energies)
    new_thetas = np.array(new_thetas)
    
    analysis = {
        'particle_name': events[0]['particles'][particle_index]['name'],
        'n_events': len(events),
        'theta_mean': np.mean(thetas),
        'theta_std': np.std(thetas),
        'theta_min': np.min(thetas),
        'new_theta_min': np.min(new_thetas),
        'theta_max': np.max(thetas),
        'new_theta_max': np.max(new_thetas),
        'energy_mean': np.mean(energies),
        'energy_std': np.std(energies),
        'energy_min': np.min(energies),
        'energy_max': np.max(energies),
        'theta_hist': np.histogram(thetas, bins=50, range=(0, 180)),
        'energy_hist': np.histogram(energies, bins=50),

    }
    
    return analysis



def transform_to_mother_frame(result):
    """
    将衰变产物变换到以母核运动方向为z轴的新坐标系
    
    参数:
    ----------
    result : dict
        breakup_kinematics函数的返回值
    
    返回:
    ----------
    dict : 包含新坐标系信息的字典
    
    字典格式:
    ----------
    {
        'mother_direction': {
            'theta_lab': float,  # 母核实验室系极角 [度]
            'phi_lab': float,    # 母核实验室系方位角 [度]
        },
        'particles': [           # 衰变粒子在新坐标系中的信息
            {
                'name': str,     # 粒子名称
                'theta_mother': float,   # 新坐标系极角 [度]，范围: 0-180
                'phi_mother': float,     # 新坐标系方位角 [度]，范围: 0-360
                'theta_lab': float,      # 实验室系极角 [度]（从原结果复制）
                'phi_lab': float,        # 实验室系方位角 [度]（从原结果复制）
            },
            ... (更多粒子)
        ],
        'angles_between': [      # 粒子两两之间的夹角 [度]
            {
                'particle1': int,        # 粒子1索引
                'particle2': int,        # 粒子2索引
                'angle': float,          # 两粒子之间的夹角 [度]
                'cos_angle': float,      # 夹角余弦值
            },
            ... (更多粒子对)
        ]
    }
    """
    
    if not result.get('valid', False):
        print("警告: 输入的衰变结果无效")
        return None
    
    # 1. 获取母核方向信息
    mother_theta_lab = result['mother_theta']
    mother_phi_lab = result['mother_phi']
    
    # 2. 创建母核方向单位向量（实验室系）
    # 母核运动方向就是新坐标系的z轴方向
    mother_theta_rad = mother_theta_lab * np.pi / 180.0
    mother_phi_rad = mother_phi_lab * np.pi / 180.0
    
    # 母核方向单位向量（实验室系）
    mother_dir_lab = np.array([
        np.sin(mother_theta_rad) * np.cos(mother_phi_rad),
        np.sin(mother_theta_rad) * np.sin(mother_phi_rad),
        np.cos(mother_theta_rad)
    ])
    
    # 3. 构建从实验室系到母核系的旋转矩阵
    # 新z轴：母核运动方向
    new_z = mother_dir_lab
    
    # 选择x轴：与原实验室系x轴的垂直分量
    lab_x = np.array([1.0, 0.0, 0.0])
    
    # 如果母核方向与x轴平行，使用y轴
    if np.abs(np.dot(new_z, lab_x)) > 0.999:
        lab_ref = np.array([0.0, 1.0, 0.0])
    else:
        lab_ref = lab_x
    
    # 新y轴：new_z × lab_ref
    new_y = np.cross(new_z, lab_ref)
    new_y = new_y / np.linalg.norm(new_y)
    
    # 新x轴：new_y × new_z
    new_x = np.cross(new_y, new_z)
    new_x = new_x / np.linalg.norm(new_x)
    
    # 旋转矩阵：从实验室系到母核系
    # 母核系坐标 = R * 实验室系坐标
    R = np.array([new_x, new_y, new_z])
    
    # 4. 变换每个衰变粒子
    particles_info = []
    
    for i, particle in enumerate(result['particles']):
        # 获取粒子动量方向单位向量（实验室系）
        p_lab = np.array([particle['px'], particle['py'], particle['pz']])
        p_mag = np.linalg.norm(p_lab)
        
        if p_mag < 1e-10:
            # 粒子静止，方向任意
            dir_lab = np.array([0.0, 0.0, 1.0])
        else:
            dir_lab = p_lab / p_mag
        
        # 变换到母核系
        dir_mother = np.dot(R, dir_lab)
        
        # 计算新坐标系中的角度
        # theta_mother: 与新z轴（母核方向）的夹角
        theta_mother = np.arccos(np.clip(dir_mother[2], -1.0, 1.0)) * 180.0 / np.pi
        
        # phi_mother: 在新x-y平面内的方位角
        if np.abs(dir_mother[0]) < 1e-10 and np.abs(dir_mother[1]) < 1e-10:
            # 沿z轴方向，phi任意，设为0
            phi_mother = 0.0
        else:
            phi_mother = np.arctan2(dir_mother[1], dir_mother[0]) * 180.0 / np.pi
        
        # 确保phi在0-360度范围内
        if phi_mother < 0:
            phi_mother += 360.0
        
        particle_info = {
            'name': particle['name'],
            'theta_mother': theta_mother,
            'phi_mother': phi_mother,
            'theta_lab': particle['theta_lab'],
            'phi_lab': particle['phi_lab'],
            'p_lab_mag': p_mag,
            'dir_mother': dir_mother.tolist(),  # 方向向量
            'index': i
        }
        
        particles_info.append(particle_info)
    
    # 5. 计算粒子两两之间的夹角
    angles_between = []
    n_particles = len(particles_info)
    
    for i in range(n_particles):
        for j in range(i+1, n_particles):
            # 获取两个粒子的方向向量
            dir_i = np.array(particles_info[i]['dir_mother'])
            dir_j = np.array(particles_info[j]['dir_mother'])
            
            # 计算夹角
            cos_angle = np.dot(dir_i, dir_j)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle = np.arccos(cos_angle) * 180.0 / np.pi
            
            angle_info = {
                'particle1': i,
                'particle1_name': particles_info[i]['name'],
                'particle2': j,
                'particle2_name': particles_info[j]['name'],
                'angle': angle,
                'cos_angle': cos_angle
            }
            
            angles_between.append(angle_info)
    
    # 6. 返回结果
    return {
        'mother_direction': {
            'theta_lab': mother_theta_lab,
            'phi_lab': mother_phi_lab,
            'dir_vector_lab': mother_dir_lab.tolist()
        },
        'coordinate_transform': {
            'new_x': new_x.tolist(),
            'new_y': new_y.tolist(),
            'new_z': new_z.tolist(),
            'rotation_matrix': R.tolist()
        },
        'particles': particles_info,
        'angles_between': angles_between,
        'n_particles': n_particles
    }




# ====================== 测试示例 ======================







# ====================== 示例衰变模式 ======================
def example_two_body_decay():
    """示例：14O* 衰变为 p + 13N"""
    print("\n示例1: 14O* → p + 13N (两体衰变)")
    
    result = breakup_kinematics(
        mass_A='14O',      # 母核
        Ek_A=14*35.0,         # 母核动能
        theta_A=0.0,      # 母核角度
        E_exciting_A=15.0,  # 母核激发能
        list_particle=['4He', '10C'],  # 衰变粒子
        list_ex=[0.0, 0.0]  # 粒子激发能
    )
    
    print_decay_result(result)
    return result


def example_four_body_decay():
    """示例：14O* 衰变为 α + α + α + 2H"""
    print("\n\n示例3: 14O* → α + α + α + 2H (四体衰变)")
    
    result = breakup_kinematics(
        mass_A='12C',
        Ek_A=12*28.0,
        theta_A=0.0,
        E_exciting_A=7.650,
        list_particle=['4He', '4He', '4He'],
        list_ex=[0.0, 0.0, 0.0]
    )
    
    print_decay_result(result)
    return result


def example_impossible_decay():
    """示例：不可能的衰变"""
    print("\n\n示例5: 不可能的衰变 (激发能不足)")
    
    result = breakup_kinematics(
        mass_A='14O',
        Ek_A=50.0,
        theta_A=0.0,
        E_exciting_A=1.0,  # 激发能太小
        list_particle=['1H', '13N'],
        list_ex=[0.0, 0.0]
    )
    
    print_decay_result(result)
    return result

def generate_decay_events_and_analyze():
    """生成多个事件并分析分布"""
    print("\n\n示例6: 生成100个衰变事件并分析分布")
    
    events = generate_multiple_events(
        mass_A='14O',
        Ek_A=14*35.0,
        theta_A=10.0,
        E_exciting_A=15.0,
        list_particle=['4He', '10C'],
        list_ex=[0.0, 0.0],
        phi_A= 10,
        n_events=10000
    )
    
    print(f"成功生成 {len(events)} 个事件")
    
    if events:
        # 分析质子分布
        proton_analysis = analyze_decay_distribution(events, particle_index=0)
        print(f"\n 4He分布分析:")
        print(f"  角度: {proton_analysis['theta_min']:.3f} - {proton_analysis['theta_max']:.3f} 度")
        print(f"  能量: {proton_analysis['energy_min']/4:.3f} - {proton_analysis['energy_max']/4:.3f} MeV/u")

        print("母核作为z轴")
        print(f"  角度: {proton_analysis['new_theta_min']:.3f} - {proton_analysis['new_theta_max']:.3f} 度")
        
        # 分析13N分布
        N13_analysis = analyze_decay_distribution(events, particle_index=1)
        print(f"\n 10C分布分析:")
        print(f"  角度: {N13_analysis['theta_min']:.3f} - {N13_analysis['theta_max']:.3f} 度")
        print(f"  能量: {N13_analysis['energy_min']/10:.3f} - {N13_analysis['energy_max']/10:.3f} MeV/u")
        print("母核作为z轴")
        print(f"  角度: {N13_analysis['new_theta_min']:.3f} - {N13_analysis['new_theta_max']:.3f} 度")

# ====================== 主函数 ======================
def main():
    """运行所有示例"""
    print("多体衰变运动学模拟器")
    print("模拟激发态原子核衰变成2-5个粒子")
    
    # 运行示例
    example_two_body_decay()
    example_four_body_decay()
    
    example_impossible_decay()
    
    # 批量生成和分析
    generate_decay_events_and_analyze()

    
    
    print("\n" + "="*80)
    print("所有示例完成")

# ====================== 高级功能 ======================




# ====================== 运行 ======================
if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)  # 批处理模式
    main()