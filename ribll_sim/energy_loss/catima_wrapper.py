"""
能损计算模块 - pycatima + ROOT TSpline3
参照 reference/brill2/range_energy_calculator 实现：
  1. 用 pycatima 预计算 1200 点射程-能量曲线
  2. 生成 TSpline3 持久化到 assets/range/{mat}_z{Z}_a{A}.root
  3. 运行时 load spline, Eval() 插值求解
"""

import os
import random
import pycatima
import ROOT


ROOT.gErrorIgnoreLevel = ROOT.kFatal


# ============================================================
#  材料配置 (与 brill2 SiliconMaterial 对应)
# ============================================================
_SPLINE_POINTS = 1200
_RANGE_EPSILON = 1e-6


def _get_material_density(mat_name):
    """获取材料密度 [g/cm³] (和 brill2 保持一致)"""
    densities = {
        'Si': 2.329,
        'C': 2.0,
    }
    return densities.get(mat_name, 1.0)


def _catima_material(mat_name):
    """获取 pycatima Material 对象"""
    z_map = {'Si': 14, 'C': 6, 'Cs': 55, 'I': 53}
    if mat_name in z_map:
        return pycatima.get_material(z_map[mat_name])
    raise ValueError("不支持的能损材料: {}".format(mat_name))


def _range_cache_path(assets_dir, mat_name, Z, A):
    """射程缓存 .root 文件路径, e.g. assets/range/Si_z6_a10.root"""
    return os.path.join(assets_dir, "range", "{mat}_z{z}_a{a}.root".format(
        mat=mat_name, z=Z, a=A))


def _assets_dir():
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "assets"
    )


# ============================================================
#  RangeEnergyCalculator — 完全参照 brill2
# ============================================================
class RangeEnergyCalculator:
    """单粒子-单材料的射程-能量插值器, 用 ROOT TSpline3"""

    def __init__(self, charge, mass, material_name, cache_path=None):
        self.charge = charge
        self.mass = mass
        self.material_name = material_name
        self.catima_mat = _catima_material(material_name)
        self.density = _get_material_density(material_name)
        self.max_energy = float(mass) * 500.0
        self.max_range = 0.0

        if cache_path is None:
            cache_path = _range_cache_path(
                _assets_dir(), material_name, charge, mass)

        self.cache_path = cache_path
        self._cache_file = None  # 保持 TFile 打开, 防止 spline 失效
        self.re_spline = None   # energy -> range
        self.er_spline = None   # range -> energy

        if not self._load_splines():
            self._build_splines()

    def _load_splines(self):
        if not os.path.exists(self.cache_path):
            return False
        f = ROOT.TFile(self.cache_path, "READ")
        if not f or f.IsZombie():
            return False
        re_spline = f.Get("re")
        er_spline = f.Get("er")
        if not re_spline or not er_spline:
            f.Close()
            return False
        self.re_spline = re_spline
        self.er_spline = er_spline
        self._cache_file = f  # 保持引用, 不关闭
        self.max_range = self.re_spline.Eval(self.max_energy)
        return True

    def _build_splines(self):
        import array
        energy = [0.0] * _SPLINE_POINTS
        rng = [0.0] * _SPLINE_POINTS

        for i in range(1, _SPLINE_POINTS):
            frac = float(i) / float(_SPLINE_POINTS - 1)
            current_energy = self.max_energy * frac * frac
            energy[i] = current_energy

            proj = pycatima.Projectile(
                float(self.mass), self.charge, 0,
                current_energy / float(self.mass)
            )
            # catima 返回 g/cm², 转为 um
            current_range_g_cm2 = pycatima.calculate(proj, self.catima_mat).range
            current_range_um = current_range_g_cm2 / self.density * 1e4

            if current_range_um <= rng[i - 1]:
                current_range_um = rng[i - 1] + _RANGE_EPSILON
            rng[i] = current_range_um

        self.max_range = rng[-1]

        arr_e = array.array('d', energy)
        arr_r = array.array('d', rng)

        self.re_spline = ROOT.TSpline3(
            "re", arr_e, arr_r, _SPLINE_POINTS,
        )
        self.er_spline = ROOT.TSpline3(
            "er", arr_r, arr_e, _SPLINE_POINTS,
        )

        cache_dir = os.path.dirname(self.cache_path)
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
        f = ROOT.TFile(self.cache_path, "RECREATE")
        self.re_spline.Write("re")
        self.er_spline.Write("er")
        f.Close()

        self._cache_file = None

    def range(self, energy):
        """能量 → 射程 [um]"""
        if energy <= 0.0:
            return 0.0
        e = min(energy, self.max_energy)
        return self.re_spline.Eval(e)

    def energy(self, rng):
        """射程 → 能量 [MeV]"""
        if rng <= 0.0:
            return 0.0
        r = min(rng, self.max_range)
        return self.er_spline.Eval(r)


# ============================================================
#  LostEnergyCalculator — 参照 brill2 LostEnergyCalculator
# ============================================================
class LostEnergyCalculator:
    """给定厚度材料的能损计算器"""

    def __init__(self, charge, mass, material_name, thickness_um, cache_path=None):
        self.calculator = RangeEnergyCalculator(charge, mass, material_name, cache_path)
        self.thickness_um = float(thickness_um)

    def residual_energy(self, incident_energy):
        """入射能量 → 出射能量"""
        if incident_energy <= 0.0:
            return 0.0
        incident_range = self.calculator.range(incident_energy)
        if incident_range <= self.thickness_um:
            return 0.0
        return max(0.0, self.calculator.energy(incident_range - self.thickness_um))

    def energy_lost(self, incident_energy):
        """入射能量 → 沉积能量"""
        if incident_energy <= 0.0:
            return 0.0
        return max(0.0, incident_energy - self.residual_energy(incident_energy))

    def incident_energy(self, residual_energy):
        """出射能量 → 入射能量 (用于靶修正)"""
        if residual_energy <= 0.0:
            return max(0.0, self.calculator.energy(self.thickness_um))
        return max(0.0, self.calculator.energy(
            self.calculator.range(residual_energy) + self.thickness_um))

    def energy_loss_from_residual(self, residual_energy):
        """出射能量 → 靶中能损"""
        if residual_energy <= 0.0:
            return max(0.0, self.incident_energy(0.0))
        return max(0.0, self.incident_energy(residual_energy) - residual_energy)


# ============================================================
#  LostEnergyCalculator 缓存池
# ============================================================
_lost_calculator_cache = {}


def _get_lost_calculator(Z, A, material_name, thickness_um):
    key = (Z, A, material_name, thickness_um)
    if key not in _lost_calculator_cache:
        _lost_calculator_cache[key] = LostEnergyCalculator(
            Z, A, material_name, thickness_um)
    return _lost_calculator_cache[key]


# ============================================================
#  对外接口 (保持和原 catima_wrapper 兼容)
# ============================================================

def calculate_energy_loss(particle_Z, particle_A, energy_MeV,
                          material_name, thickness_um):
    """
    计算粒子穿过材料后的剩余动能 (CSDA, TSpline3 插值)

    参数:
        particle_Z: 原子序数
        particle_A: 质量数
        energy_MeV: 入射动能 [MeV]
        material_name: 材料名 ('C', 'Si')
        thickness_um: 厚度 [um]

    返回:
        float: 剩余动能 [MeV], 阻停则返回 0.0
    """
    if energy_MeV <= 0.0 or thickness_um <= 0.0:
        return energy_MeV

    calc = _get_lost_calculator(particle_Z, particle_A, material_name, thickness_um)
    return calc.residual_energy(energy_MeV)


def calculate_range(particle_Z, particle_A, energy_MeV, material_name):
    """计算粒子在材料中的射程 [um]"""
    if energy_MeV <= 0.0:
        return 0.0
    return _get_lost_calculator(
        particle_Z, particle_A, material_name, 0.0).calculator.range(energy_MeV)


def will_particle_stop(particle_Z, particle_A, energy_MeV,
                       material_name, thickness_um):
    """判断粒子是否阻停在材料中"""
    return calculate_range(
        particle_Z, particle_A, energy_MeV, material_name) <= thickness_um


def half_target_correction(particle_Z, particle_A, energy_MeV_in_detector,
                           target_material='C', target_thickness_um=100.0):
    """
    半靶厚能损修正: 从探测器测到的能量反推反应能量
    假设粒子在靶中心反应, 出靶损失半个靶厚能量
    """
    half_thickness = target_thickness_um / 2.0
    calc = RangeEnergyCalculator(particle_Z, particle_A, target_material)

    if energy_MeV_in_detector <= 0.0:
        return max(0.0, calc.energy(half_thickness))

    # 迭代: E_reaction = E_detected + EnergyLostInTarget(E_reaction)
    e_guess = energy_MeV_in_detector + 0.5
    for _ in range(10):
        r = calc.range(e_guess)
        if r <= half_thickness:
            loss = e_guess
        else:
            loss = e_guess - calc.energy(r - half_thickness)
        new_guess = energy_MeV_in_detector + loss
        if abs(new_guess - e_guess) < 1e-6:
            break
        e_guess = new_guess

    return e_guess