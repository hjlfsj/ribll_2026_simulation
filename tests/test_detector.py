"""
探测器几何模块测试
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from ribll_sim.detector.geometry import Target, DSSD, SSD, DetectorArray


def test_target_sampling():
    """测试靶内均匀采样"""
    print("测试1: 靶内反应点采样")
    target = Target(radius=15.0, thickness=100.0)

    n_inside = 0
    n_total = 1000
    for _ in range(n_total):
        pt = target.sample_reaction_point()
        x, y, z = pt['true']
        if target.is_inside(x, y, z):
            n_inside += 1

    print(f"  靶内比例: {n_inside}/{n_total}")
    assert n_inside == n_total, "有点在靶外!"
    print("  通过\n")


def test_target_ppac_smear():
    """测试PPAC位置模糊"""
    print("测试2: PPAC位置模糊")
    target = Target(radius=15.0, thickness=100.0)

    for _ in range(10):
        pt = target.sample_reaction_point(ppac_sigma_x=1.0, ppac_sigma_y=1.0)
        _, _, z_obs = pt['observed']
        assert z_obs == 0.0, f"z_obs应为0, 实际{z_obs}"

    print("  z_obs始终为0 (正确)")
    print("  通过\n")


def test_dssd_grid():
    """测试DSSD网格"""
    print("测试3: DSSD网格定位")
    dssd = DSSD('test', z_position=100.0, thickness=1000, grid_n=64)

    cx, cy, i, j = dssd.get_grid_center(0.0, 0.0)
    print(f"  击中(0,0) -> 网格中心({cx:.3f}, {cy:.3f}), 索引({i},{j})")
    assert 31 <= i <= 32 and 31 <= j <= 32, f"中心网格索引异常: ({i},{j})"

    cx, cy, i, j = dssd.get_grid_center(31.5, 31.5)
    print(f"  击中(31.5,31.5) -> 网格中心({cx:.3f}, {cy:.3f}), 索引({i},{j})")
    assert i == 63 and j == 63, f"边角网格索引异常: ({i},{j})"

    print("  通过\n")


def test_detector_hit():
    """测试粒子击中探测器判定"""
    print("测试4: 粒子击中判定")
    array = DetectorArray().build_default()

    hits = array.trace_particle(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    print(f"  正入射(0,0,1): 击中{len(hits)}个探测器")
    assert len(hits) == 5, f"应击中5个, 实际{len(hits)}"

    hits = array.trace_particle(0.0, 0.0, 0.0, 0.0, 1.0, 1.0)
    print(f"  斜入射(0,1,1): 击中{len(hits)}个探测器")

    hits = array.trace_particle(0.0, 0.0, 0.0, 0.0, 10.0, 1.0)
    print(f"  大角度入射(0,10,1): 击中{len(hits)}个探测器")
    print("  通过\n")


def test_detector_array_order():
    """测试探测器z排序"""
    print("测试5: 探测器z排序")
    array = DetectorArray().build_default()

    z_list = [d.z for d in array.detectors]
    print(f"  z位置: {z_list}")
    assert z_list == sorted(z_list), "探测器未按z排序"
    assert len(array.detectors) == 5
    assert array.detectors[0].name == 'DSSD1'
    assert array.detectors[-1].name == 'SSD'
    print("  通过\n")


if __name__ == '__main__':
    print("=" * 60)
    print("探测器几何模块测试")
    print("=" * 60 + "\n")

    test_target_sampling()
    test_target_ppac_smear()
    test_dssd_grid()
    test_detector_hit()
    test_detector_array_order()

    print("=" * 60)
    print("全部测试完成")
    print("=" * 60)