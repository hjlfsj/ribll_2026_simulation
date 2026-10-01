"""
能损计算测试
计算粒子在硅(Si)和碳(C)中的射程、阻停判断及剩余能量
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import time

from ribll_sim.energy_loss.catima_wrapper import (
    calculate_energy_loss,
    calculate_range,
    will_particle_stop,
)


def print_header(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def run_case(name, Z, A, energies, material, thickness_um):
    """对一组能量计算并打印结果"""
    print(f"\n{'─' * 70}")
    print(f"  {name} (Z={Z}, A={A}) → {material} {thickness_um}um")
    print(f"  {'能量 [MeV]':<14} {'射程 [um]':<14} {'阻停':<8} {'剩余能量 [MeV]'}")
    print(f"  {'─' * 14} {'─' * 14} {'─' * 8} {'─' * 20}")

    for energy in energies:
        rng = calculate_range(Z, A, energy, material)
        stops = will_particle_stop(Z, A, energy, material, thickness_um)
        residual = calculate_energy_loss(Z, A, energy, material, thickness_um)
        print(f"  {energy:<14.2f} {rng:<14.1f} {str(stops):<8} {residual:<.2f}")


def test_dssd_layers():
    """各 DSSD/SSD 层厚度下的能损"""
    print_header("探测器各层 Si 能损 (按默认探测器阵列厚度)")

    layers = [
        ("DSSD1", 68),
        ("DSSD2-4", 1000),
        ("SSD", 1500),
    ]

    particles = [
        ("10C", 6, 10, [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]),
        ("4He", 2, 4,  [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]),
        ("6He", 2, 6,  [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]),
        ("14O", 8, 14, [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]),
        ("6Li", 3, 6,  [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]),
    ]

    for layer_name, thickness in layers:
        print(f"\n  >>> {layer_name} ({thickness}um Si) <<<")
        for name, Z, A, energies in particles:
            run_case(name, Z, A, energies, "Si", thickness)


def test_target_C():
    """靶 (100um C) 中的能损"""
    print_header("靶中能损 (C, 100um)")

    energies = [1, 5, 10, 20, 50, 100, 150, 200, 300, 500]

    particles = [
        ("10C", 6, 10),
        ("4He", 2, 4),
        ("6He", 2, 6),
        ("14O", 8, 14),
        ("6Li", 3, 6),
    ]

    for name, Z, A in particles:
        run_case(name, Z, A, energies, "C", 100)


def test_threshold():
    """二分法求临界穿透能量"""
    print_header("临界穿透能量 (恰好穿透指定厚度)")

    cases = [
        ("10C", 6, 10, "Si", 68,   "DSSD1"),
        ("10C", 6, 10, "Si", 1000, "DSSD2-4"),
        ("10C", 6, 10, "Si", 1500, "SSD"),
        ("4He", 2, 4,  "Si", 68,   "DSSD1"),
        ("4He", 2, 4,  "Si", 1000, "DSSD2-4"),
        ("4He", 2, 4,  "Si", 1500, "SSD"),
        ("6He", 2, 6,  "Si", 68,   "DSSD1"),
        ("6He", 2, 6,  "Si", 1000, "DSSD2-4"),
        ("6He", 2, 6,  "Si", 1500, "SSD"),
        ("14O", 8, 14, "Si", 68,   "DSSD1"),
        ("14O", 8, 14, "Si", 1000, "DSSD2-4"),
        ("14O", 8, 14, "Si", 1500, "SSD"),
        ("6Li", 3, 6,  "Si", 68,   "DSSD1"),
        ("6Li", 3, 6,  "Si", 1000, "DSSD2-4"),
        ("6Li", 3, 6,  "Si", 1500, "SSD"),
    ]

    print(f"\n  {'粒子':<8} {'探测器':<12} {'厚度 [um]':<12} {'临界能量 [MeV]':<18} {'该能量射程 [um]'}")
    print(f"  {'─'*8} {'─'*12} {'─'*12} {'─'*18} {'─'*20}")

    for name, Z, A, mat, thick, det_name in cases:
        lo, hi = 0.01, 1000.0
        for _ in range(35):
            mid = (lo + hi) / 2
            if will_particle_stop(Z, A, mid, mat, thick):
                lo = mid
            else:
                hi = mid
        E_th = hi
        rng = calculate_range(Z, A, E_th, mat)
        print(f"  {name:<8} {det_name:<12} {thick:<12} {E_th:<18.4f} {rng:<.1f}")


def test_performance():
    """批量计算性能"""
    print_header("批量计算性能")

    Z, A, mat, thick = 6, 10, "Si", 1000
    energies = [i * 10.0 for i in range(1, 1001)]

    t0 = time.perf_counter()
    results = [calculate_energy_loss(Z, A, e, mat, thick) for e in energies]
    t1 = time.perf_counter()

    n_stopped = sum(1 for r in results if r == 0.0)
    print(f"\n  10C in Si 1000um, 能量范围 10-10000 MeV, 共 1000 次")
    print(f"  耗时: {(t1-t0)*1000:.1f} ms  ({1000/(t1-t0):.0f} 次/秒)")
    print(f"  阻停: {n_stopped}/1000,  穿透: {1000-n_stopped}/1000")


if __name__ == '__main__':
    print("=" * 70)
    print("  能损计算模块测试")
    print("=" * 70)

    test_dssd_layers()
    test_target_C()
    test_threshold()
    test_performance()

    print()
    print("=" * 70)
    print("  全部完成")
    print("=" * 70)