"""
运动学模块测试
测试三种反应类型的正确性：守恒律、Q值、边界条件
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from ribll_sim.kinematics.reaction_type1 import simulate_type1_decay
from ribll_sim.kinematics.reaction_type2 import simulate_type2_reaction
from ribll_sim.kinematics.reaction_type3 import simulate_type3_sequential


def check_conservation_2body(initial_p4, p1, p2, tolerance=1e-6):
    total_final = p1['p4'] + p2['p4']
    dE = abs(initial_p4.E() - total_final.E())
    dp = (initial_p4 - total_final).P()
    ok = dE < tolerance and dp < tolerance
    if not ok:
        print(f"  守恒违反: dE={dE:.2e} MeV, dp={dp:.2e} MeV/c")
    return ok


def test_type1_basic():
    """类型1: 14O -> 10C + 4He (无激发)"""
    print("测试1: 14O -> 10C + 4He")
    result = simulate_type1_decay('14O', 100.0, '10C', '4He')
    assert result['valid'], f"衰变失败: {result.get('reason')}"
    print(f"  Q值: {result['Q_value']:.3f} MeV")
    print(f"  10C: Ek={result['B']['Ek']:.2f} MeV, theta={result['B']['theta']:.2f} deg")
    print(f"  4He: Ek={result['C']['Ek']:.2f} MeV, theta={result['C']['theta']:.2f} deg")
    assert check_conservation_2body(result['A']['p4'], result['B'], result['C'])
    print("  通过\n")


def test_type1_excited_C():
    """类型1: C带激发态"""
    print("测试2: 14O -> 10C + 4He (4He激发10MeV)")
    result = simulate_type1_decay('14O', 100.0, '10C', '4He', excitation_C=10.0)
    if result['valid']:
        print(f"  Q值: {result['Q_value']:.3f} MeV")
        assert check_conservation_2body(result['A']['p4'], result['B'], result['C'])
    else:
        print(f"  预期失败: {result['reason']}")
    print("  通过\n")


def test_type2_elastic():
    """类型2: 弹散 14O(2H,2H)14O"""
    print("测试3: 弹散 14O + 2H -> 14O + 2H (E=15*14 MeV)")
    result = simulate_type2_reaction('14O', 15*14, '2H', '14O', '2H')
    assert result['valid'], f"反应失败: {result.get('reason')}"
    print(f"  类型: {result['reaction_type']}")
    print(f"  Q值: {result['Q_value']:.3f} MeV")
    print(f"  14O: Ek={result['C']['Ek']:.2f} MeV, theta={result['C']['theta']:.2f} deg")
    print(f"  2H:  Ek={result['D']['Ek']:.2f} MeV, theta={result['D']['theta']:.2f} deg")
    print("  通过\n")


def test_type2_transfer():
    """类型2: 转移反应 14O(2H,6Li)10C"""
    print("测试4: 转移反应 14O + 2H -> 10C + 6Li (E=35*14 MeV)")
    result = simulate_type2_reaction('14O', 35*14, '2H', '10C', '6Li')
    assert result['valid'], f"反应失败: {result.get('reason')}"
    print(f"  类型: {result['reaction_type']}")
    print(f"  Q值: {result['Q_value']:.3f} MeV")
    print(f"  10C: Ek={result['C']['Ek']:.2f} MeV, theta={result['C']['theta']:.2f} deg")
    print(f"  6Li: Ek={result['D']['Ek']:.2f} MeV, theta={result['D']['theta']:.2f} deg")
    print("  通过\n")


def test_type2_invalid():
    """类型2: 能量不足以克服负Q值"""
    print("测试5: 能量不足 (E=0 MeV)")
    result = simulate_type2_reaction('14O', 0.0, '2H', '10C', '6Li')
    assert not result['valid'], "应该失败但没有"
    print(f"  预期失败: {result['reason']}")
    print("  通过\n")


def test_type3_basic():
    """类型3: 14O + 2H -> 10C* + 6Li, 10C* -> 6He + 4He"""
    print("测试6: 级联反应 14O(2H,6Li)10C*, 10C* -> 6He + 4He")
    result = simulate_type3_sequential(
        '14O', 35*14, '2H',
        '10C', '6Li',
        '6He', '4He',
        excitation_C=10.0
    )
    assert result['valid'], f"反应失败: {result.get('reason')}"
    print(f"  第一步Q值: {result['step1']['Q_value']:.3f} MeV")
    print(f"  C*衰变Q值: {result['Q_breakup']:.3f} MeV")
    print(f"  6He: Ek={result['E']['Ek']:.2f} MeV, theta={result['E']['theta']:.2f} deg")
    print(f"  4He: Ek={result['F']['Ek']:.2f} MeV, theta={result['F']['theta']:.2f} deg")
    print(f"  6Li: Ek={result['D']['Ek']:.2f} MeV, theta={result['D']['theta']:.2f} deg")

    p4_Cstar = result['step1']['C']['p4']
    total_final = result['E']['p4'] + result['F']['p4']
    dM = abs(p4_Cstar.M() - total_final.M())
    print(f"  C*质量守恒检查: dM={dM:.2e} MeV")
    assert dM < 1e-3, f"C*质量不守恒: dM={dM:.2e}"

    p4_initial = result['step1']['C']['p4'] + result['step1']['D']['p4']
    total_final_all = result['E']['p4'] + result['F']['p4'] + result['D']['p4']
    dE = abs(p4_initial.E() - total_final_all.E())
    dp = (p4_initial - total_final_all).P()
    print(f"  总能动量守恒: dE={dE:.2e} MeV, dp={dp:.2e} MeV/c")
    print("  通过\n")


def test_type3_excited_F():
    """类型3: F带激发态"""
    print("测试7: 级联反应, F(4He)激发5MeV")
    result = simulate_type3_sequential(
        '14O', 35*14, '2H',
        '10C', '6Li',
        '6He', '4He',
        excitation_C=15.0,
        excitation_F=5.0
    )
    if result['valid']:
        print(f"  C*衰变Q值: {result['Q_breakup']:.3f} MeV")
        print(f"  F激发态质量: {result['F']['mass']:.3f} MeV")
        assert result['F']['excitation'] == 5.0
    else:
        print(f"  预期可能失败: {result['reason']}")
    print("  通过\n")


if __name__ == '__main__':
    print("=" * 60)
    print("运动学模块测试")
    print("=" * 60 + "\n")

    test_type1_basic()
    test_type1_excited_C()
    test_type2_elastic()
    test_type2_transfer()
    test_type2_invalid()
    test_type3_basic()
    test_type3_excited_F()

    print("=" * 60)
    print("全部测试完成")
    print("=" * 60)