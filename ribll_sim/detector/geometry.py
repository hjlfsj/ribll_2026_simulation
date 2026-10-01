"""
探测器几何模块
可复用的探测器类，通过config.toml配置参数

探测器类型：
- Target: 圆柱形靶
- DSSD: 双面硅条探测器 (有网格)
- SSD: 硅探测器 (无网格)
- CsI: 碘化铯探测器 (预留)
"""

import random
from math import sqrt


class Target:
    """圆柱形反应靶，对称轴为z轴，几何中心在(0,0,0)"""

    def __init__(self, radius=15.0, thickness=100.0):
        """
        参数:
        ----------
        radius : float
            靶半径 [mm], 默认15mm (直径30mm)
        thickness : float
            靶厚度 [um], 默认100um
        """
        self.radius = radius
        self.thickness_um = thickness
        self.thickness_mm = thickness * 1e-3
        self.material = 'C'  # 默认碳靶

    def sample_reaction_point(self, ppac_sigma_x=0.0, ppac_sigma_y=0.0):
        """
        在靶体积内均匀随机生成反应点，并叠加PPAC位置模糊

        参数:
        ----------
        ppac_sigma_x : float
            PPAC x方向位置分辨 [mm] (sigma), 默认0
        ppac_sigma_y : float
            PPAC y方向位置分辨 [mm] (sigma), 默认0

        返回:
        ----------
        tuple : (x, y, z) [mm], z始终在 [-thickness_mm/2, thickness_mm/2] 均匀分布
        """
        r = self.radius * sqrt(random.random())
        phi = random.random() * 2 * 3.141592653589793
        x_true = r * sqrt(1 - random.random()) * (1 if random.random() < 0.5 else -1)
        x_true = r * (2 * random.random() - 1)
        y_true = r * (2 * random.random() - 1)
        while x_true**2 + y_true**2 > self.radius**2:
            x_true = r * (2 * random.random() - 1)
            y_true = r * (2 * random.random() - 1)

        z_true = (random.random() - 0.5) * self.thickness_mm

        x_obs = x_true + random.gauss(0, ppac_sigma_x) if ppac_sigma_x > 0 else x_true
        y_obs = y_true + random.gauss(0, ppac_sigma_y) if ppac_sigma_y > 0 else y_true
        z_obs = 0.0  # 实验上z位置不可知

        return {
            'true': (x_true, y_true, z_true),
            'observed': (x_obs, y_obs, z_obs)
        }

    def is_inside(self, x, y, z):
        """判断点是否在靶内"""
        half = self.thickness_mm / 2.0
        return (x**2 + y**2 <= self.radius**2) and (-half <= z <= half)


class Detector:
    """探测器基类"""

    def __init__(self, name, z_position, thickness, width=64.0, height=64.0):
        """
        参数:
        ----------
        name : str
            探测器名称
        z_position : float
            探测器平面z坐标 [mm]
        thickness : float
            探测器厚度 [um]
        width : float
            探测器x方向宽度 [mm], 默认64
        height : float
            探测器y方向高度 [mm], 默认64
        """
        self.name = name
        self.z = z_position
        self.thickness_um = thickness
        self.thickness_mm = thickness * 1e-3
        self.width = width
        self.height = height
        self.half_w = width / 2.0
        self.half_h = height / 2.0

    def can_hit(self, x, y, z, px, py, pz):
        """
        判断从(x,y,z)出发、动量为(px,py,pz)的粒子能否打到探测器平面

        返回:
        ----------
        tuple : (can_hit: bool, hit_x: float, hit_y: float)
        """
        if pz <= 0:
            return (False, 0.0, 0.0)

        t = (self.z - z) / pz
        hit_x = x + px * t
        hit_y = y + py * t

        if (abs(hit_x) <= self.half_w and abs(hit_y) <= self.half_h):
            return (True, hit_x, hit_y)

        return (False, hit_x, hit_y)

    def path_length_in_detector(self, px, py, pz):
        """粒子穿过探测器的路径长度 [mm]，假设垂直入射近似"""
        if pz <= 0:
            return float('inf')
        total_p = sqrt(px**2 + py**2 + pz**2)
        cos_theta = pz / total_p
        if cos_theta <= 0:
            return float('inf')
        return self.thickness_mm / cos_theta


class SiliconDetector(Detector):
    """硅探测器基类"""

    def __init__(self, name, z_position, thickness, width=64.0, height=64.0):
        super().__init__(name, z_position, thickness, width, height)
        self.material = 'Si'


class DSSD(SiliconDetector):
    """双面硅条探测器，具有网格结构"""

    def __init__(self, name, z_position, thickness, grid_n=64, width=64.0, height=64.0):
        """
        参数:
        ----------
        grid_n : int
            网格数 (64或32)
        """
        super().__init__(name, z_position, thickness, width, height)
        self.grid_n = grid_n
        self.grid_size = width / grid_n  # mm per grid

    def get_grid_center(self, hit_x, hit_y):
        """
        粒子击中位置对应的网格中心坐标

        参数:
        ----------
        hit_x, hit_y : float
            击中位置 [mm]

        返回:
        ----------
        tuple : (grid_center_x, grid_center_y, grid_i, grid_j)
        """
        i = int((hit_x + self.half_w) / self.grid_size)
        j = int((hit_y + self.half_h) / self.grid_size)

        i = max(0, min(self.grid_n - 1, i))
        j = max(0, min(self.grid_n - 1, j))

        cx = -self.half_w + (i + 0.5) * self.grid_size
        cy = -self.half_h + (j + 0.5) * self.grid_size

        return (cx, cy, i, j)

    def get_grid_index(self, hit_x, hit_y):
        """获取网格索引"""
        i = int((hit_x + self.half_w) / self.grid_size)
        j = int((hit_y + self.half_h) / self.grid_size)
        i = max(0, min(self.grid_n - 1, i))
        j = max(0, min(self.grid_n - 1, j))
        return (i, j)


class SSD(SiliconDetector):
    """普通硅探测器，无网格"""
    pass


class CsIDetector(Detector):
    """CsI探测器 (预留)"""

    def __init__(self, name, z_position, thickness, width=64.0, height=64.0):
        super().__init__(name, z_position, thickness, width, height)
        self.material = 'CsI'


class DetectorArray:
    """探测器阵列管理器"""

    def __init__(self, target=None):
        self.target = target if target else Target()
        self.detectors = []

    def add_detector(self, detector):
        self.detectors.append(detector)
        self.detectors.sort(key=lambda d: d.z)

    def build_default(self):
        """构建默认的RIBLL2026探测器配置"""
        self.target = Target(radius=15.0, thickness=100.0)
        self.detectors = [
            DSSD('DSSD1', z_position=100.0, thickness=68, grid_n=32),
            DSSD('DSSD2', z_position=110.0, thickness=1000, grid_n=64),
            DSSD('DSSD3', z_position=120.0, thickness=1000, grid_n=32),
            DSSD('DSSD4', z_position=130.0, thickness=1000, grid_n=32),
            SSD('SSD', z_position=140.0, thickness=1500),
            CsIDetector('CsI', z_position=200.0, thickness=40000, width=80.0, height=80.0),
        ]
        return self

    def trace_particle(self, origin_x, origin_y, origin_z, px, py, pz):
        """
        追踪从靶点出发的粒子，返回依次击中的探测器列表

        返回:
        ----------
        list[dict] : 每个元素包含 detector, hit_x, hit_y
        """
        hits = []
        for det in self.detectors:
            can, hx, hy = det.can_hit(origin_x, origin_y, origin_z, px, py, pz)
            if can:
                hits.append({
                    'detector': det,
                    'hit_x': hx,
                    'hit_y': hy,
                })
        return hits