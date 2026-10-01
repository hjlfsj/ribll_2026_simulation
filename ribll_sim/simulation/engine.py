"""
RIBLL2026 模拟主引擎
负责：加载配置 -> 生成反应 -> 追踪粒子 -> 计算能损 -> 记录事件
"""

import os
import sys
import time
import random
from math import sqrt
from datetime import datetime

import ROOT
import numpy as np

from ..data.nuclear_data import get_particle_mass, parse_particle
from ..kinematics import simulate_type1_decay, simulate_type2_reaction, simulate_type3_sequential
from ..detector.geometry import Target, DSSD, SSD, DetectorArray
from ..energy_loss.catima_wrapper import calculate_energy_loss, will_particle_stop


def parse_config(config):
    """从config.toml字典中提取参数"""
    reaction = config['reaction']
    det_cfg = config.get('detectors', {})
    sim_cfg = config.get('simulation', {})

    return {
        'reaction_type': reaction['type'],
        'particle_A': reaction.get('particle_A', '14O'),
        'particle_B': reaction.get('particle_B', '2H'),
        'particle_C': reaction.get('particle_C', '10C'),
        'particle_D': reaction.get('particle_D', '6Li'),
        'particle_E': reaction.get('particle_E', '6He'),
        'particle_F': reaction.get('particle_F', '4He'),
        'E_beam': reaction.get('E_beam', 35.0 * 14.0),
        'excitation_C': reaction.get('excitation_C', 0.0),
        'excitation_F': reaction.get('excitation_F', 0.0),
        'cms_theta': reaction.get('cms_theta', 180.0),
        'n_events': sim_cfg.get('n_events', 10000),
        'target_thickness': det_cfg.get('target_thickness_um', 100.0),
        'ppac_sigma': det_cfg.get('ppac_sigma_xy', 0.0),
        'si_resolution': det_cfg.get('si_resolution_fwhm', 0.01),
    }


def build_detector_array(config):
    """根据config构建探测器阵列"""
    det_cfg = config.get('detectors', {})
    target_thickness = det_cfg.get('target_thickness_um', 100.0)

    array = DetectorArray(target=Target(radius=15.0, thickness=target_thickness))

    detector_list = det_cfg.get('detector_list', None)
    if detector_list is None:
        return array.build_default()

    for d in detector_list:
        dtype = d['type']
        name = d['name']
        z = d['z']
        thickness = d['thickness']
        width = d.get('width', 64.0)
        height = d.get('height', 64.0)

        if dtype == 'DSSD':
            grid_n = d.get('grid_n', 64)
            array.add_detector(DSSD(name, z, thickness, grid_n, width, height))
        elif dtype == 'SSD':
            array.add_detector(SSD(name, z, thickness, width, height))

    return array


def run_simulation(config):
    """运行模拟主循环"""
    params = parse_config(config)
    array = build_detector_array(config)

    print("=" * 60)
    print("RIBLL2026 模拟")
    print("=" * 60)
    print(f"反应类型: {params['reaction_type']}")
    print(f"粒子: {params['particle_A']} + {params['particle_B']}")
    print(f"束流能量: {params['E_beam']:.1f} MeV")
    print(f"事件数: {params['n_events']}")
    print(f"靶厚度: {params['target_thickness']} um")
    print(f"Si分辨率(FWHM): {params['si_resolution']*100:.1f}%")
    print(f"探测器: {[d.name for d in array.detectors]}")
    print()

    n_valid = 0
    start_time = time.time()

    results_theory = []
    results_detected = []

    for i in range(params['n_events']):
        if i % max(1, params['n_events'] // 10) == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            print(f"  进度: {i}/{params['n_events']} ({100*i/params['n_events']:.0f}%), "
                  f"速率: {rate:.0f} evt/s, 有效: {n_valid}")

        kin_result = generate_event(params)
        if kin_result is None or not kin_result.get('valid'):
            continue

        particles = extract_particles(kin_result, params['reaction_type'])

        theory_event = record_theory_event(particles, kin_result, params['reaction_type'])
        results_theory.append(theory_event)

        detected = process_detector_response(particles, array, params)
        if detected is not None:
            n_valid += 1
            results_detected.append(detected)

    total_time = time.time() - start_time
    print(f"\n完成! 总时间: {total_time:.1f}s, "
          f"理论事件: {len(results_theory)}, 探测事件: {n_valid}")

    return {
        'theory': results_theory,
        'detected': results_detected,
        'params': params,
        'n_total': params['n_events'],
        'n_valid': n_valid,
        'time': total_time
    }


def generate_event(params):
    """根据反应类型生成单个事件"""
    rt = params['reaction_type']

    if rt == 1:
        return simulate_type1_decay(
            params['particle_A'], params['E_beam'],
            params['particle_B'], params['particle_C'],
            excitation_C=params['excitation_C']
        )
    elif rt == 2:
        return simulate_type2_reaction(
            params['particle_A'], params['E_beam'], params['particle_B'],
            params['particle_C'], params['particle_D'],
            excitation_C=params['excitation_C'],
            cms_theta=params['cms_theta']
        )
    elif rt == 3:
        return simulate_type3_sequential(
            params['particle_A'], params['E_beam'], params['particle_B'],
            params['particle_C'], params['particle_D'],
            params['particle_E'], params['particle_F'],
            excitation_C=params['excitation_C'],
            excitation_F=params['excitation_F'],
            cms_theta=params['cms_theta']
        )
    return None


def extract_particles(kin_result, reaction_type):
    """从运动学结果提取粒子列表"""
    particles = []
    if reaction_type == 1:
        particles = [kin_result['B'], kin_result['C']]
    elif reaction_type == 2:
        particles = [kin_result['C'], kin_result['D']]
    elif reaction_type == 3:
        particles = [kin_result['E'], kin_result['F'], kin_result['D']]
    return particles


def record_theory_event(particles, kin_result, reaction_type):
    """记录理论事件（精确运动学值）"""
    event = {
        'reaction_type': reaction_type,
        'Q_value': kin_result.get('Q_value', kin_result.get('step1', {}).get('Q_value', 0)),
        'n_particles': len(particles),
    }
    for i, p in enumerate(particles):
        event[f'p{i}_name'] = p['name']
        event[f'p{i}_Ek'] = p['Ek']
        event[f'p{i}_theta'] = p['theta']
        event[f'p{i}_phi'] = p['phi']
        event[f'p{i}_p'] = p['p']
        event[f'p{i}_mass'] = p['mass']
    return event


def process_detector_response(particles, array, params):
    """
    模拟探测器响应：
    1. 靶中能损 -> 出靶
    2. 依次穿过各探测器(Si+CsI)，计算能损
    3. 粒子必须被探测到(Si或CsI阻停)
    4. 对测量值加分辨模糊(Si=1%FWHM, CsI=3%FWHM)
    """
    target = array.target
    si_resolution = params['si_resolution']
    csi_resolution = 0.03
    ppac_sigma = params['ppac_sigma']

    pt = target.sample_reaction_point(ppac_sigma, ppac_sigma)
    ox, oy, oz = pt['true']

    detected_particles = []
    all_stopped = True

    for p in particles:
        name = p['name']
        p4 = p['p4']
        px, py, pz = p4.Px(), p4.Py(), p4.Pz()
        Ek = p['Ek']
        Z, A = parse_particle(name)

        if pz <= 0:
            all_stopped = False
            break

        Ek_after_target = calculate_energy_loss(Z, A, Ek, 'C', target.thickness_um / 2.0)
        if Ek_after_target <= 0:
            all_stopped = False
            break

        total_detected_E = 0.0
        stopped_in_csi = False
        stopped_detector_name = None
        remaining_E = Ek_after_target

        for det in array.detectors:
            can_hit, hx, hy = det.can_hit(ox, oy, oz, px, py, pz)
            if not can_hit:
                continue

            if isinstance(det, DSSD) or isinstance(det, SSD):
                E_after_det = calculate_energy_loss(Z, A, remaining_E, 'Si', det.thickness_um)
                energy_deposited = remaining_E - E_after_det
                total_detected_E += energy_deposited

                if E_after_det <= 0:
                    stopped_detector_name = det.name
                    remaining_E = 0
                    break
                else:
                    remaining_E = E_after_det

            elif hasattr(det, 'material') and det.material == 'CsI':
                stopped_in_csi = True
                stopped_detector_name = det.name
                total_detected_E += remaining_E
                remaining_E = 0
                break

        if stopped_detector_name is None:
            all_stopped = False
            break

        detected_particles.append({
            'name': name,
            'Ek_true': Ek,
            'Ek_detected': total_detected_E,
            'theta_true': p['theta'],
            'stopped_in': stopped_detector_name,
            'stopped_in_csi': stopped_in_csi,
        })

    if not all_stopped:
        return None

    event = {
        'reaction_point': (ox, oy, oz),
        'n_particles': len(detected_particles),
    }
    for i, dp in enumerate(detected_particles):
        res = csi_resolution if dp['stopped_in_csi'] else si_resolution
        sigma = res * dp['Ek_detected'] / 2.355
        Ek_smeared = dp['Ek_detected'] + random.gauss(0, sigma) if sigma > 0 else dp['Ek_detected']

        event[f'det_p{i}_name'] = dp['name']
        event[f'det_p{i}_Ek'] = max(0, Ek_smeared)
        event[f'det_p{i}_Ek_true'] = dp['Ek_detected']
        event[f'det_p{i}_theta_true'] = dp['theta_true']
        event[f'det_p{i}_stopped_in'] = dp['stopped_in']

    return event


def save_to_root(results, output_dir, timestamp):
    """保存模拟结果到ROOT文件"""
    os.makedirs(output_dir, exist_ok=True)

    root_path = os.path.join(output_dir, f'sim_{timestamp}.root')
    f = ROOT.TFile(root_path, 'RECREATE')

    theory_data = results['theory']
    n_theory = len(theory_data)

    if n_theory > 0:
        keys = list(theory_data[0].keys())
        tree_theory = ROOT.TTree('theory', 'Theoretical Events')

        arrays = {}
        for key in keys:
            val = theory_data[0][key]
            if isinstance(val, float):
                arrays[key] = np.zeros(1, dtype=np.float64)
                tree_theory.Branch(key, arrays[key], f'{key}/D')
            elif isinstance(val, int):
                arrays[key] = np.zeros(1, dtype=np.int32)
                tree_theory.Branch(key, arrays[key], f'{key}/I')

        for evt in theory_data:
            for key in keys:
                arrays[key][0] = evt[key]
            tree_theory.Fill()

        tree_theory.Write()

    detected_data = results['detected']
    n_detected = len(detected_data)

    if n_detected > 0:
        keys_d = list(detected_data[0].keys())
        keys_d = [k for k in keys_d if k != 'reaction_point']
        tree_det = ROOT.TTree('detected', 'Detected Events')

        arrays_d = {}
        for key in keys_d:
            val = detected_data[0][key]
            if isinstance(val, float):
                arrays_d[key] = np.zeros(1, dtype=np.float64)
                tree_det.Branch(key, arrays_d[key], f'{key}/D')
            elif isinstance(val, int):
                arrays_d[key] = np.zeros(1, dtype=np.int32)
                tree_det.Branch(key, arrays_d[key], f'{key}/I')
            elif isinstance(val, str):
                continue

        for evt in detected_data:
            for key in keys_d:
                if key in arrays_d:
                    arrays_d[key][0] = evt[key]
            tree_det.Fill()

        tree_det.Write()

    info = ROOT.TNamed('info', f'total:{results["n_total"]} '
                               f'valid:{results["n_valid"]} '
                               f'time:{results["time"]:.1f}s')
    info.Write()

    f.Close()
    print(f"\n结果已保存至: {root_path}")
    print(f"  理论事件: {n_theory}")
    print(f"  探测事件: {n_detected}")

    return root_path