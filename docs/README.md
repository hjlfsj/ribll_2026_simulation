# RIBLL2026 模拟框架文档

## 1. 项目概述

本项目是一个基于 ROOT + Python 的核物理模拟框架，用于模拟 RIBLL2026 实验中的核反应过程。实现从运动学计算、探测器响应模拟、能量损失修正到数据可视化的完整模拟链。

- **工作目录**: `/home/ribll2026/ribll2026_www/github_code/ribll_2026_simulation`
- **输出目录**: `/data/disk1/ribll2026_www_data/simulation`
- **语言**: Python 3 (依赖 ROOT、pycatima、tomli、PyQt6、pyvista、matplotlib)

## 2. 目录结构

```
ribll_2026_simulation/
├── simulation.py              # 主入口脚本
├── config.toml                # 模拟参数配置文件
├── assets/                    # 静态数据资源
│   ├── nuclear_data.py        # 核素质量数据（与ribll_sim/data中一致）
│   ├── range/                 # 预计算的射程-能量 TSpline3 缓存
│   │   ├── Si_z6_a10.root     # 10C 在 Si 中的射程曲线
│   │   ├── Si_z3_a6.root      # 6Li 在 Si 中的射程曲线
│   │   ├── C_z6_a10.root      # 10C 在 C 中的射程曲线
│   │   └── ...                # 其他粒子-材料组合
│   └── *.txt                  # 验证用理论计算数据
├── docs/
│   └── README.md              # 本文档
├── reference/                 # 参考实现
│   └── pyroot/
│       ├── standard_tools/    # ROOT版参考代码
│       │   ├── generate_kinematics.py
│       │   └── nuclear_data.py
│       └── test_module/       # 参考测试代码
│           ├── code/
│           │   └── test_type_abcd.py
│           ├── data_root/     # 验证用 .root 文件
│           └── data_pdf/      # 验证用 .pdf 曲线
└── ribll_sim/                 # 核心模拟包
    ├── __init__.py
    ├── data/
    │   ├── __init__.py
    │   └── nuclear_data.py    # 核素质量数据库
    ├── kinematics/
    │   ├── __init__.py
    │   ├── reaction_type1.py  # A → B + C 衰变
    │   ├── reaction_type2.py  # A + B → C + D 反应
    │   ├── reaction_type3.py  # A + B → C* + D, C* → E + F
    │   └── reconstruction.py  # 激发能重建（理论与实验）
    ├── detector/
    │   ├── __init__.py
    │   └── geometry.py        # 探测器几何与响应
    ├── energy_loss/
    │   ├── __init__.py
    │   └── catima_wrapper.py  # pycatima 封装（TSpline3 缓存）
    ├── simulation/
    │   ├── __init__.py
    │   └── engine.py          # 批量模拟引擎
    └── gui/
        ├── __init__.py
        └── app.py             # GUI 模块（PyQt6 + pyvista + matplotlib）
```

## 3. 物理模型

### 3.1 三种反应类型

所有运动学计算均基于相对论四动量 (`ROOT.TLorentzVector`)。Type 1 使用 `ROOT.TGenPhaseSpace` 进行均匀相空间生成；Type 2 和 Type 3 第一步使用质心系角度采样的确定性方法。

#### 类型 1: A → B + C（衰变）

母核 A 以动能 $E_k$ 沿 z 轴飞行，在飞行中衰变为 B 和 C。采用各向同性相空间 (`TGenPhaseSpace`)。

| 参数 | 说明 |
|---|---|
| `particle_A` | 母核符号，如 `14O` |
| `E_beam` | 母核动能 [MeV] |
| `particle_B` | 衰变产物 B（始终基态） |
| `particle_C` | 衰变产物 C（可设激发态） |
| `excitation_C` | C 的激发能 [MeV] → 物理上等价于母核 A 激发 |

**激发能规则**: Type 1 仅有 `excitation_C` 有效。设置 `excitation_C > 0` 相当于母核 A 处于激发态（$m_A = m_{A,gs} + E_x$），衰变产物 C 获得额外质量。

**质心系限制**: Type 1 **不使用** `cms_theta` 参数，衰变始终各向同性。

#### 类型 2: A + B → C + D（转移反应）

入射粒子 A 以动能 $E_k$ 撞击静止靶核 B，生成 C 和 D。使用质心系角度采样方法（非 `TGenPhaseSpace`）。

| 类别 | 条件 |
|---|---|
| 弹性散射 | C = A, D = B, excitation_C = 0 |
| 非弹性散射 | C = A, D = B, excitation_C > 0 |
| 转移反应 | C ≠ A 或 D ≠ B |

| 参数 | 说明 |
|---|---|
| `particle_A` | 入射粒子 |
| `particle_B` | 靶核 |
| `particle_C` | 出射粒子 C（可设激发态） |
| `particle_D` | 出射粒子 D（始终基态） |
| `excitation_C` | C 的激发能 [MeV] |
| `cms_theta` | C 在 A+B 质心系中的最大极角 [度] |

**激发能规则**: Type 2 仅有 `excitation_C` 有效（产物 C 激发）。`particle_D` 始终处于基态。

**质心系限制**: C 在 A+B 质心系中的极角 $\theta_{cm}$ 在 $(0, \texttt{cms\_theta})$ 范围内均匀抽样，方位角 $\phi_{cm}$ 在 $(0, 2\pi)$ 范围内均匀抽样。默认 `cms_theta = 180`（无限制，全覆盖）。

**运动学计算流程**:
1. 计算质心系总能量 $\sqrt{s} = M(P_A + P_B)$ 和质心系动量 $p_{cm} = \frac{\sqrt{[s-(m_C+m_D)^2][s-(m_C-m_D)^2]}}{2\sqrt{s}}$
2. 在 $(0, \texttt{cms\_theta})$ 内均匀抽样 $\theta_{cm}$，在 $(0, 2\pi)$ 内均匀抽样 $\phi_{cm}$
3. 在质心系中构建 C 和 D 的四动量（C 沿采样方向，D 反向）
4. Boost 回实验室系得到最终四动量
5. 返回 `cms_theta_C` 字段记录 C 的实际质心系角度

**判断动能是否足够**:
- 若 Q 值 < 0，需要阈能 $E_{threshold} = -Q \cdot (m_A + m_B) / m_B$，$E_{beam}$ 必须大于此值

#### 类型 3: A + B → C\* + D, C\* → E + F（级联反应）

先进行类型 2 反应生成激发态 C\* 和 D，随后 C\* 衰变为 E 和 F。

| 步骤 | 方法 | 角度分布 |
|---|---|---|
| 第一步: A+B → C\*+D | 质心系角度采样（同 Type 2） | C\* 在 $(0, \texttt{cms\_theta})$ 内均匀 |
| 第二步: C\* → E+F | `TGenPhaseSpace` | C\* 静止系中各向同性 |

| 参数 | 说明 |
|---|---|
| `particle_E` | C\* 衰变产物 E（始终基态） |
| `particle_F` | C\* 衰变产物 F（可设激发态） |
| `excitation_C` | C 的激发能 [MeV]（用于第一步，生成 C\*） |
| `excitation_F` | F 的激发能 [MeV]（用于第二步，F 可激发） |
| `cms_theta` | 第一步 C\* 在 A+B 质心系中的最大极角 [度] |

**激发能规则**: Type 3 中 `excitation_C` 用于中间态 C\*（第一步），`excitation_F` 用于衰变产物 F（第二步）。`particle_D` 和 `particle_E` 始终基态。

**质心系限制**: 仅第一步使用 `cms_theta` 限制（同 Type 2）。第二步 C\* → E+F 衰变始终各向同性（`TGenPhaseSpace`）。

**坐标变换细节**: E、F 的方向先在 C\* 静止系中（以 C\* 方向为 z 轴）通过 `TGenPhaseSpace` 生成，再变换到实验室系。

### 3.2 产物粒子信息结构

每个产物粒子返回以下字段：

```python
{
    'name': str,           # 粒子符号 (如 '10C')
    'mass_gs': float,      # 基态质量 [MeV/c²]
    'excitation': float,   # 激发能 [MeV]
    'mass': float,         # 实际质量 (gs + excitation) [MeV/c²]
    'p4': TLorentzVector,  # ROOT 四动量
    'Ek': float,           # 动能 [MeV]
    'theta': float,        # 极角 [度]
    'phi': float,          # 方位角 [度]
    'p': float,            # 动量大小 [MeV/c]
}
```

类型 2 结果额外包含:
```python
{
    'Q_value': float,          # 反应 Q 值 [MeV]
    'reaction_type': str,      # 'elastic' / 'inelastic' / 'transfer'
    'cms_theta_C': float,      # C 在质心系中的实际极角 [度]
}
```

类型 3 结果额外包含:
```python
{
    'step1': dict,             # 第一步 (Type 2) 的完整结果
    'cms_theta_C': float,      # C* 在第一步质心系中的实际极角 [度]
    'Q_breakup': float,        # C* 衰变 Q 值 [MeV]
}
```

### 3.3 激发能重建

#### 理论重建 (`reconstruct_excitation`)

使用运动学结果的精确四动量（不含探测器效应）重建激发能：

| 反应类型 | 重建目标 | 公式 |
|---|---|---|
| 类型 1 | A | $E_x(A) = M_{reco}(B+C) - m_{A,gs}$ |
| 类型 2 | C | $E_x(C) = M(P_A+P_B-P_D) - m_{C,gs}$（缺失质量法） |
| 类型 3 | C | $E_x(C) = M(E+F) - m_{C,gs}$（不变质量法） |

#### 实验重建 (`reconstruct_excitation_experimental`)

使用探测器测量值（包含能量分辨模糊和角度重建误差）重建激发能。该函数使用靶校正后的实验能量 `Eexp`（而非探测沉积能量）和实验测量角度 $\theta$ 构建四动量。

| 反应类型 | 所用数据 | 粒子名称来源 |
|---|---|---|
| 类型 1 | det_p0 (B) + det_p1 (C) 的 Eexp, θ | `params['particle_A/B/C']` |
| 类型 2 | det_p1 (D) 的 Eexp, θ（缺失质量法） | `params['particle_C/D']` |
| 类型 3 | det_p0 (E) + det_p1 (F) 的 Eexp, θ（不变质量法） | `params['particle_C/E/F']` |

**关键**: 实验重建使用 `Eexp`（靶校正后能量 = 探测沉积能量总和 + 靶中能损），确保反映粒子在反应发生时的真实动能。

## 4. 能量损失模型

### 4.1 TSpline3 缓存架构 (参照 brill2)

能量损失计算参照 `reference/brill2/range_energy_calculator` 实现，使用预计算的 ROOT `TSpline3` 曲线进行插值，避免实时调用 `pycatima` 导致的性能问题。

**核心类**:

| 类 | 职责 |
|---|---|
| `RangeEnergyCalculator(charge, mass, material_name)` | 单粒子-单材料的射程↔能量双向插值器 |
| `LostEnergyCalculator(charge, mass, material_name, thickness_um)` | 给定厚度的能损计算器 |

**工作流程**:

1. **首次运行**: 对每个 (粒子 Z,A, 材料) 组合，调用 `pycatima` 计算 1200 个 quadratically-spaced 能量点的射程，生成两个 `TSpline3` 曲线:
   - `re`: energy → range（能量→射程）
   - `er`: range → energy（射程→能量）
   - 持久化到 `assets/range/{mat}_z{Z}_a{A}.root`

2. **后续运行**: 直接从 `.root` 文件加载 `TSpline3`，通过 `Eval()` 插值求解

3. **能损计算**: 
   - `residual_energy(E_in)` = `er(range(E_in) - thickness)`：入射→出射
   - `energy_lost(E_in)` = `E_in - residual_energy(E_in)`：入射→沉积
   - `incident_energy(E_out)` = `er(range(E_out) + thickness)`：出射→入射（靶修正用）

**缓存位置**: `assets/range/` 目录，文件命名格式 `{材料}_{z}{Z}_a{A}.root`，如:
- `Si_z6_a10.root` — 10C 在 Si 中的射程曲线
- `C_z3_a6.root` — 6Li 在 C (靶) 中的射程曲线

**支持的材料**:

| 材料名 | pycatima Z | 密度 [g/cm³] |
|---|---|---|
| C | 6 | 2.0 |
| Si | 14 | 2.329 |
| Cs | 55 | — |
| I | 53 | — |

### 4.2 对外接口

| 函数 | 说明 |
|---|---|
| `calculate_energy_loss(Z, A, E, material, thickness_um)` | 穿过材料后剩余能量 [MeV]。阻停则返回 0.0 |
| `calculate_range(Z, A, E, material)` | 粒子在材料中射程 [μm] |
| `will_particle_stop(Z, A, E, material, thickness_um)` | 判断粒子是否阻停在材料中 |
| `half_target_correction(Z, A, E_det, material, thickness_um)` | 半靶厚能损修正：从探测器测能反推反应能量（迭代法） |

## 5. 探测器模型

### 5.1 几何定义

- **坐标原点 (0,0,0)** 为靶的几何中心
- **z 轴** 为束流方向（也是圆柱靶对称轴）
- 入射粒子从 z < 0 方向射入

### 5.2 探测器类型

| 类型 | 尺寸 | 材料 | 说明 |
|---|---|---|---|
| `Target` | 直径 30mm, 厚默认 100μm | C（碳靶） | 圆柱形，反应点均匀随机分布 |
| `DSSD` | 64×64 mm² | Si | 双面硅条，带网格 (32×32 或 64×64) |
| `SSD` | 64×64 mm² | Si | 普通硅探测器，无网格 |
| `CsIDetector` | 80×80 mm² | CsI | 碘化铯探测器，用于阻停穿透 Si 的高能粒子 |

### 5.3 默认探测器阵列

| 名称 | 类型 | z 位置 | 厚度 | 网格 | 能量分辨 (FWHM) |
|---|---|---|---|---|---|
| DSSD1 | DSSD | 100 mm | 68 μm | 32×32 | 1% |
| DSSD2 | DSSD | 110 mm | 1000 μm | 64×64 | 1% |
| DSSD3 | DSSD | 120 mm | 1000 μm | 32×32 | 1% |
| DSSD4 | DSSD | 130 mm | 1000 μm | 32×32 | 1% |
| SSD | SSD | 140 mm | 1500 μm | — | 1% |
| CsI | CsI | 200 mm | 40000 μm (40mm) | — | 3% |

### 5.4 模拟探测流程（GUI 模式详细版）

以下为 GUI 主窗口中 `_process_detector` 的完整模拟流程：

1. **反应点采样**: 在靶体积内均匀随机生成反应点 $(x,y,z)$，PPAC 位置模糊叠加在观测 $(x,y)$ 上

2. **靶中能损（半靶厚模型）**: 假设粒子在靶中心 (z≈0) 反应，出靶损失半个靶厚的能量。
   - `target_loss = Ek - Ek_after_target`
   - 半靶厚 = `target_thickness / 2`

3. **逐层探测器追踪**: 对每个探测器依次判断粒子能否命中（`can_hit`），若命中则计算能损。

4. **Si 探测器层 (DSSD1-4, SSD) 处理**:
   - 计算理论能损: `deposit = remaining_E - residual_energy(remaining_E, 'Si', thickness)`
   - 叠加能量模糊: `deposit_smeared = deposit + Gauss(0, sigma)`，其中 `sigma = si_resolution * deposit / 2.355`
   - 若 `residual_energy ≤ 0`，粒子在该层阻停。若该层是 DSSD，角度重建使用网格中心
   - 若粒子未阻停，继续下一层

5. **CsI 探测器层处理**:
   - 剩余能量全部沉积: `deposit = remaining_E`
   - 叠加能量模糊: `sigma = 3% * deposit / 2.355`
   - 粒子在 CsI 中阻停

6. **粒子未阻停处理**:
   - 若粒子穿透所有探测器（Si+CsI），标记为 "离开探测器阵列"，事件无效

7. **每层信息记录**（用于 GUI 分层能损对比展示）:
   ```python
   {
       'name': str,           # 层名 (如 'DSSD2')
       'material': str,       # 'Si' 或 'CsI'
       'thickness_um': float,
       'E_theory': float,     # 理论能损 [MeV]
       'E_detected': float,   # 探测能损 [MeV]（叠加模糊后）
       'remaining': float,    # 剩余能量 [MeV]
       'stopped': bool,       # 是否在该层阻停
   }
   ```

8. **靶校正实验能量**:
   - `Eexp = total_detected_E + target_loss`
   - 即探测器测得的总沉积能量 + 靶中损失的能量

9. **实验激发能重建**: 调用 `reconstruct_excitation_experimental`，使用 `Eexp` 和测量角度 $\theta$

### 5.5 角度重建

- **DSSD（粒子阻停在 DSSD 中）**: 使用粒子击中的网格中心 $(cx, cy)$ 和 PPAC 观测位置 $(ox_{obs}, oy_{obs})$ 反算:
  $$\theta_{reco} = \arctan\left(\frac{\sqrt{(cx - ox_{obs})^2 + (cy - oy_{obs})^2}}{z_{det}}\right)$$

- **SSD / CsI（粒子阻停在非 DSSD 中）**: 使用实际击中位置 $(hx, hy)$ 和 PPAC 观测位置反算

- **粒子穿透未阻停**: 使用理论角度 $\theta_{theory}$

## 6. 配置文件

`config.toml` 是模拟的唯一参数来源，支持热修改无需重编译。

```toml
# RIBLL2026 模拟配置文件

# ====== 激发能设置规则 ======
# Type 1 (A→B+C):     只有 A 能设置激发能 → excitation_C 用于产物 C
#                       (物理上 A*→B+C, C 质量=基态质量+激发能)
# Type 2 (A+B→C+D):   只有 C 能设置激发能 → excitation_C 用于产物 C
# Type 3 (A+B→C*+D→E+F): 只有 C 和 F 能设置激发能 → excitation_C 用于 C*,
#                         excitation_F 用于衰变产物 F
#                         D 和 E 始终基态 (不存在 excitation_D)
#
# ====== 质心系角度限制 ======
# cms_theta: 产物 C 在 A+B 质心系中的最大极角 [度]
#            第一步反应中, C 的质心系角度在 (0, cms_theta) 中均匀抽样
#            默认 180 度 (无限制, 全覆盖)
#            Type 1 不使用此参数 (衰变始终各向同性)
#            Type 3 第二步 C*→E+F 衰变始终各向同性 (TGenPhaseSpace)

[reaction]
# 反应类型: 1=A→B+C, 2=A+B→C+D, 3=A+B→C*+D, C*→E+F
type = 2

# 入射粒子
particle_A = "14O"
E_beam = 490.0              # MeV (35 MeV/u * 14)

# 靶核
particle_B = "2H"

# 出射粒子 (类型2和3使用)
particle_C = "10C"
particle_D = "6Li"

# C* 衰变产物 (仅类型3使用)
particle_E = "6He"
particle_F = "4He"

# 激发能 [MeV]
# Type 1: 只有 excitation_C 有效 (A 激发)
# Type 2: 只有 excitation_C 有效 (C 激发)
# Type 3: excitation_C (C*激发) + excitation_F (F 激发) 有效
excitation_C = 3.3
excitation_F = 0.0

# 质心系角度限制 [度] (仅 Type 2,3 第一步使用)
cms_theta = 20.0

[detectors]
# 靶参数
target_thickness_um = 100.0

# PPAC 位置模糊: x-y 平面 sigma [mm]
ppac_sigma_xy = 0.0

# Si 探测器能量分辨 FWHM (比例)
si_resolution_fwhm = 0.01

# 自定义探测器列表（可选，不配置则使用默认阵列）
[[detectors.detector_list]]
name = "DSSD1"
type = "DSSD"
z = 100.0
thickness = 68
grid_n = 32

[[detectors.detector_list]]
name = "DSSD2"
type = "DSSD"
z = 110.0
thickness = 1000
grid_n = 64

[[detectors.detector_list]]
name = "DSSD3"
type = "DSSD"
z = 120.0
thickness = 1000
grid_n = 32

[[detectors.detector_list]]
name = "DSSD4"
type = "DSSD"
z = 130.0
thickness = 1000
grid_n = 32

[[detectors.detector_list]]
name = "SSD"
type = "SSD"
z = 140.0
thickness = 1500

[simulation]
n_events = 10000

[output]
data_dir = "/data/disk1/ribll2026_www_data/simulation"
```

### 支持的粒子符号

`1H`, `2H`, `3H`, `4H`, `3He`, `4He`, `5He`, `6He`, `8He`, `6Li`, `7Li`, `8Li`, `9Li`, `11Li`, `7Be`, `9Be`, `10Be`, `11Be`, `10B`, `11B`, `12B`, `10C`, `11C`, `12C`, `13C`, `14C`, `15C`, `13N`, `14N`, `15N`, `14O`, `15O`, `16O`, `17O`, `18O`

## 7. 主程序用法

### 7.1 启动 GUI 模式（默认）

```bash
cd /home/ribll2026/ribll2026_www/github_code/ribll_2026_simulation
python3 simulation.py
```

GUI 包含三个窗口：

#### 窗口 1: 主窗口（3D 可视化 + 事件信息）

| 组件 | 内容 |
|---|---|
| **右上角 Quit 按钮** | 关闭整个程序 |
| **控制面板** | `Batch Mode` 复选框、`N` 事件数输入框、`Sim` 按钮、`Save .root` 按钮 |
| **状态标签** | 显示模拟状态 |
| **左侧信息面板** (占 33%) | 事件详情（见下方） |
| **右侧 3D 视图** (占 67%) | 探测器几何 (靶+探测器+DSSD 网格) + 粒子径迹 + z 轴参考线 |

**左侧信息面板显示内容**:

```
=== Event #N ===
Type: R  —  反应描述
cms_theta: X°  (C in CMS: 0–X°)
Beam: 14O @ 490.0 MeV
Target: 2H
Q-value: -11.941 MeV
E_x(10C): 3.300 MeV
Reaction point: x=X.XXX mm, y=X.XXX mm, z=X.X um

产物粒子1: Ek=XXX.XX MeV, theta=X.XXdeg, phi=X.XXdeg
产物粒子2: Ek=XXX.XX MeV, theta=X.XXdeg, phi=X.XXdeg
C CMS angle: XX.XX°

--- 10C (Z=6, A=10)  Ek0=424.28 MeV ---  ← 粒子1
  Layer            ThE-loss  DetE-loss     Remain     状态
  ---------------------------------------------------
  C靶/2              4.542         --    419.742     穿过     ← 半靶厚
  DSSD1               6.460      6.436    413.282     穿过     ← 低饱和度蓝色
  DSSD2             106.516    106.327    306.766     穿过     ← 低饱和度绿色
  DSSD3             153.524    153.603    153.242     穿过     ← 低饱和度紫色
  DSSD4             153.242    153.135      0.000     阻停     ← 低饱和度红色
  ---------------------------------------------------
   Total(Si+CsI)    419.742    419.501  测量theta=2.77°
           E_exp    424.043  (靶校正后粒子能量)
  阻停在: DSSD4

--- 6Li (Z=3, A=6)  Ek0=53.77 MeV ---    ← 粒子2
  ...

========================================================
E_x(10C) 设置值:   3.300 MeV      ← config.toml 中 excitation_C
E_x(10C) 理论重建: 3.300 MeV      ← 精确四动量重建
E_x(10C) 实验重建: 3.285 MeV      ← 探测器测量值重建（含模糊效应）
```

**3D 视图特性**:

| 元素 | 描述 |
|---|---|
| 探测器几何 | 靶（黄色圆柱）、各 Si 探测器（低饱和度不同颜色平板）、DSSD 网格线、CsI 探测器 |
| z 轴参考线 | 灰色虚线，线宽 3，范围 z ∈ (-10mm, 200mm) |
| 粒子径迹 | 彩色实线，从反应点出发到阻停位置（或透射到远处） |
| 粒子标注 | 在 z=50mm 处标注粒子名称 |
| 阻停标记 | 粒子阻停位置绘制彩色小球 |

**探测器层颜色**（低饱和度）:

| 层 | 颜色 | 色值 |
|---|---|---|
| C靶/2 | 低饱和度金 | `#C8A96E` |
| DSSD1 | 低饱和度蓝 | `#9DC3E6` |
| DSSD2 | 低饱和度绿 | `#A9D18E` |
| DSSD3 | 低饱和度紫 | `#B4A0D4` |
| DSSD4 | 低饱和度红 | `#D4A0A0` |
| SSD | 低饱和度青 | `#A0D4D0` |
| CsI | 低饱和度橙 | `#D4C0A0` |

#### 窗口 2: 理论分析窗口（Theory Analysis）

- 2×2 布局: 各粒子理论 E-θ 二维直方图 + 理论激发能谱
- y 轴标签: $E_{particle}$ [MeV]
- x 轴标签: $\theta_{particle}$ [deg]
- 激发能谱可交互 rebin（调整 bins, min, max）
- 数据来源: 运动学精确四动量，无探测器效应

#### 窗口 3: 实验分析窗口（Detected Analysis）

- 2×2 布局: 各粒子实验 E-θ 二维直方图 + 实验激发能谱
- y 轴标签: $E_{particle}^{exp}$ [MeV]（靶校正后实验能量）
- x 轴标签: $\theta_{particle}$ [deg]
- 激发能谱可交互 rebin
- 数据来源: 探测器测量值（含能量分辨模糊和角度重建误差）
- 使用 `det_p{i}_Eexp`（靶校正后能量）和 `det_p{i}_theta`（测量角度）
- 使用 `Ex_exp`（实验重建激发能）

### 7.2 批量模拟模式

```bash
python3 simulation.py -n 10000
```

跳过 GUI，直接模拟 10000 个事件，输出 ROOT 文件。

### 7.3 命令行参数

| 参数 | 说明 | 示例 |
|---|---|---|
| `-n N` | 模拟 N 个事件（非 GUI 模式） | `-n 50000` |
| `-c PATH` | 指定配置文件路径 | `-c my_config.toml` |

## 8. 测试程序用法

### 8.1 测试运动学模块

在 Python 交互环境中直接测试各运动学函数：

```python
import sys
sys.path.insert(0, '/home/ribll2026/ribll2026_www/github_code/ribll_2026_simulation')

from ribll_sim.kinematics import (
    simulate_type1_decay,
    simulate_type2_reaction,
    simulate_type3_sequential,
)
from ribll_sim.kinematics.reconstruction import (
    reconstruct_excitation,
    reconstruct_excitation_experimental,
)

# 测试类型 1: A -> B + C 衰变
result = simulate_type1_decay(
    particle_A='14O', E_A_kin=490.0,
    particle_B='6Li', particle_C='10C',
    excitation_C=3.3
)
print("类型1 valid:", result['valid'])
print(f"  B: Ek={result['B']['Ek']:.1f} MeV, theta={result['B']['theta']:.1f}°")
print(f"  C: Ek={result['C']['Ek']:.1f} MeV, theta={result['C']['theta']:.1f}°")

# 测试类型 2: 转移反应（含质心系角度限制）
result = simulate_type2_reaction(
    particle_A='14O', E_A_kin=490.0, particle_B='2H',
    particle_C='10C', particle_D='6Li',
    excitation_C=3.3,
    cms_theta=20.0,          # 质心系角度限制
)
print("类型2 valid:", result['valid'])
print(f"  Q_value: {result['Q_value']:.3f} MeV")
print(f"  反应类别: {result['reaction_type']}")
print(f"  C CMS angle: {result['cms_theta_C']:.2f}°")

# 测试类型 3: 级联反应（含质心系角度限制）
result = simulate_type3_sequential(
    particle_A='14O', E_A_kin=490.0, particle_B='2H',
    particle_C='10C', particle_D='6Li',
    particle_E='6He', particle_F='4He',
    excitation_C=3.3, excitation_F=2.0,
    cms_theta=20.0,          # 第一步质心系角度限制
)
print("类型3 valid:", result['valid'])
if result.get('step1'):
    print(f"  Step1 Q: {result['step1']['Q_value']:.3f} MeV")
    print(f"  C* CMS angle: {result['cms_theta_C']:.2f}°")
```

### 8.2 测试激发能重建（理论与实验对比）

```python
from ribll_sim.kinematics import simulate_type2_reaction
from ribll_sim.kinematics.reconstruction import (
    reconstruct_excitation,
    reconstruct_excitation_experimental,
)

params = {
    'reaction_type': 2,
    'particle_A': '14O', 'E_beam': 490.0,
    'particle_B': '2H', 'particle_C': '10C', 'particle_D': '6Li',
}

# 生成事件并获取运动学结果
result = simulate_type2_reaction(
    params['particle_A'], params['E_beam'], params['particle_B'],
    params['particle_C'], params['particle_D'],
    excitation_C=3.3, cms_theta=20.0,
)

# 理论重建（使用精确四动量）
ex_theory = reconstruct_excitation(result, params)
print(f"真实 E_x=3.3 MeV, 理论重建 E_x={ex_theory:.3f} MeV")

# 模拟实验重建（需要构造 detected_data 模拟探测器响应后）
# 实际使用见 gui/app.py 中的 _process_detector + reconstruct_excitation_experimental
```

### 8.3 测试探测器几何

```python
from ribll_sim.detector.geometry import Target, DSSD, DetectorArray

target = Target(radius=15.0, thickness=100.0)
pt = target.sample_reaction_point(ppac_sigma_x=1.0, ppac_sigma_y=1.0)
print(f"真实位置: {pt['true']}")
print(f"PPAC 观测: {pt['observed']}")

# 构建默认阵列（含 CsI）
array = DetectorArray(target=target).build_default()
print(f"探测器列表: {[d.name for d in array.detectors]}")
# 输出: ['DSSD1', 'DSSD2', 'DSSD3', 'DSSD4', 'SSD', 'CsI']

# 追踪粒子
hits = array.trace_particle(
    origin_x=pt['true'][0], origin_y=pt['true'][1],
    origin_z=pt['true'][2], px=10.0, py=5.0, pz=100.0,
)
for h in hits:
    d = h['detector']
    print(f"命中 {d.name} @ ({h['hit_x']:.1f}, {h['hit_y']:.1f}) @ z={d.z}")

# 测试 DSSD 网格定位
if hits and isinstance(hits[0]['detector'], DSSD):
    dssd = hits[0]['detector']
    cx, cy, gi, gj = dssd.get_grid_center(hits[0]['hit_x'], hits[0]['hit_y'])
    print(f"网格中心: ({cx:.1f}, {cy:.1f}), 网格索引: ({gi}, {gj})")
```

### 8.4 测试能量损失

```python
from ribll_sim.energy_loss.catima_wrapper import (
    calculate_energy_loss, calculate_range, half_target_correction,
)

# 10C 离子穿过 50μm 碳靶（半靶厚）
E_remaining = calculate_energy_loss(
    particle_Z=6, particle_A=10, energy_MeV=424.28,
    material_name='C', thickness_um=50.0,
)
print(f"10C 入射 424.28 MeV → 半靶厚后剩余: {E_remaining:.2f} MeV")

# 计算射程
r = calculate_range(particle_Z=6, particle_A=10, energy_MeV=424.0, material_name='Si')
print(f"10C 在 Si 中射程: {r:.1f} μm")

# 半靶厚修正（探测器测能 45 MeV → 靶前反应能）
E_corrected = half_target_correction(
    particle_Z=3, particle_A=6, energy_MeV_in_detector=45.0,
    target_material='C', target_thickness_um=100.0,
)
print(f"探测器测能 45 MeV → 靶前反应能: {E_corrected:.2f} MeV")

# 首次运行后检查缓存文件
# ls assets/range/
# 应看到 Si_z6_a10.root、C_z6_a10.root 等文件
```

### 8.5 参考测试代码

参考代码位于 `reference/pyroot/test_module/code/test_type_abcd.py`，可使用以下命令运行：

```bash
cd /home/ribll2026/ribll2026_www/github_code/ribll_2026_simulation
python3 reference/pyroot/test_module/code/test_type_abcd.py
```

## 9. 模块 API 速查

### 9.1 `ribll_sim.data.nuclear_data`

| 函数 | 说明 |
|---|---|
| `get_particle_mass(particle)` | 获取粒子基态质量 [MeV/c²] |
| `parse_particle(particle)` | 解析粒子符号为 `(Z, A)` |
| `get_mass_data()` | 返回所有核素质量数据库字典 |

**质量计算公式**: $m = \text{mass\_excess\_keV}[1] + A \times 931.494\ \text{MeV}$

### 9.2 `ribll_sim.kinematics`

| 函数 | 说明 |
|---|---|
| `simulate_type1_decay(A, E, B, C, excitation_C=0.0)` | 类型 1 衰变模拟（各向同性 TGenPhaseSpace） |
| `simulate_type2_reaction(A, E, B, C, D, excitation_C=0.0, cms_theta=180.0)` | 类型 2 反应模拟（质心系角度采样） |
| `simulate_type3_sequential(A, E, B, C, D, E_p, F_p, excitation_C=0.0, excitation_F=0.0, cms_theta=180.0)` | 类型 3 级联反应模拟 |
| `reconstruct_excitation(kin_result, params)` | 自动选择方法重建理论激发能 |
| `reconstruct_excitation_experimental(detected_data, kin_result, params)` | 使用探测器数据重建实验激发能（含分辨模糊） |

### 9.3 `ribll_sim.detector.geometry`

| 类 | 说明 |
|---|---|
| `Target(radius=15.0, thickness=100.0)` | 圆柱靶，`sample_reaction_point(sigma_x, sigma_y)` 采样反应点 |
| `DSSD(name, z, thickness, grid_n=64, width=64, height=64)` | 双面硅条，带 `get_grid_center(hx, hy)` 和 `get_grid_index(hx, hy)` |
| `SSD(name, z, thickness, width=64, height=64)` | 普通硅探测器 |
| `CsIDetector(name, z, thickness, width=80, height=80)` | CsI 探测器（material='CsI'，3% 能量分辨） |
| `DetectorArray(target=None)` | 探测器阵列，`build_default()` 构建默认配置，`trace_particle()` 追踪粒子 |

### 9.4 `ribll_sim.energy_loss.catima_wrapper`

| 类/函数 | 说明 |
|---|---|
| `RangeEnergyCalculator(charge, mass, material_name, cache_path)` | 射程↔能量双向插值器（TSpline3） |
| `LostEnergyCalculator(charge, mass, material_name, thickness_um)` | 给定厚度的能损计算器 |
| `calculate_energy_loss(Z, A, E, material, thickness_um)` | 穿过材料后剩余能量 |
| `calculate_range(Z, A, E, material)` | 粒子在材料中射程 [μm] |
| `will_particle_stop(Z, A, E, material, thickness_um)` | 判断粒子是否阻停 |
| `half_target_correction(Z, A, E_det, material, thickness_um)` | 半靶厚能损修正（迭代法） |

### 9.5 `ribll_sim.simulation.engine`

| 函数 | 说明 |
|---|---|
| `run_simulation(config)` | 批量模拟，返回 `{'theory': [...], 'detected': [...], 'params': ..., 'n_total': ..., 'n_valid': ..., 'time': ...}` |
| `save_to_root(results, output_dir, timestamp)` | 保存结果为 ROOT 文件（TTree） |
| `parse_config(config)` | 从 config.toml 字典提取参数 |
| `generate_event(params)` | 根据反应类型生成单个运动学事件 |
| `process_detector_response(particles, array, params)` | 模拟探测器响应 |

### 9.6 `ribll_sim.gui.app`

| 类/函数 | 说明 |
|---|---|
| `MainWindow(config)` | 主窗口：3D 可视化 + 事件信息面板 + 控制按钮 |
| `AnalysisWindow(title, params, product_names, is_detected)` | 分析窗口：E-θ 二维图 + 激发能谱（理论/实验） |
| `RIBLLApp(config)` | 应用控制器：创建三个窗口并管理模拟逻辑 |
| `parse_config_simple(config)` | 简化版配置解析（GUI 用） |
| `get_product_names(params)` | 根据反应类型返回产物粒子列表 |
| `get_excited_particle_name(params)` | 返回被激发粒子名称（用于显示标签） |

## 10. 输出格式

每次批量模拟在输出目录下创建一个时间戳命名的子目录：

```
/data/disk1/ribll2026_www_data/simulation/
└── 20260930_143052/
    ├── simulation.root     # TTree: theory + detected
    └── config.toml         # 本次模拟参数快照
```

ROOT 文件包含两个 TTree:

**`theory` TTree** — 理论事件（精确运动学，无探测器效应）:

| Branch | 类型 | 说明 |
|---|---|---|
| `reaction_type` | I | 反应类型 |
| `Q_value` | D | Q 值 [MeV] |
| `n_particles` | I | 产物粒子数 |
| `p{i}_name` | — | 粒子名称 (string, 不存入 TTree) |
| `p{i}_Ek` | D | 动能 [MeV] |
| `p{i}_theta` | D | 极角 [度] |
| `p{i}_phi` | D | 方位角 [度] |
| `p{i}_p` | D | 动量大小 [MeV/c] |
| `p{i}_mass` | D | 质量 [MeV/c²] |

**`detected` TTree** — 探测事件（含能量分辨模糊）:

| Branch | 类型 | 说明 |
|---|---|---|
| `n_particles` | I | 产物粒子数 |
| `det_p{i}_name` | — | 粒子名称 (string, 不存入 TTree) |
| `det_p{i}_Ek` | D | 探测能量（叠加模糊后）[MeV] |
| `det_p{i}_Ek_true` | D | 真实沉积能量（未模糊）[MeV] |
| `det_p{i}_theta_true` | D | 理论角度 [度] |
| `det_p{i}_stopped_in` | — | 阻停探测器名 (string, 不存入 TTree) |

**`info` TNamed**: 元信息（总事件数、有效事件数、模拟耗时）

## 11. GUI 保存格式（Save .root）

在 GUI 中点击 `Save .root` 按钮，保存为 TH2D/TH1D 格式：

ROOT 文件中的对象:
- `theory_h2_p0`, `theory_h2_p1`, ... — 各粒子理论 E-θ 二维谱
- `theory_h_ex` — 理论激发能谱（精确重建）
- `detected_h2_p0`, `detected_h2_p1`, ... — 各粒子实验 E-θ 二维谱（使用 `Eexp` 和测量 $\theta$）
- `detected_h_ex` — 实验激发能谱（探测器数据重建）

## 12. 依赖环境

| 依赖 | 版本要求 | 用途 |
|---|---|---|
| Python | ≥ 3.8 | 运行环境 |
| ROOT | ≥ 6.24 | 四动量计算（TLorentzVector）、TSpline3、TTree、TGenPhaseSpace |
| pycatima | ≥ 1.3 | 能量损失基础计算（仅用于生成 TSpline3 缓存，运行时不再调用） |
| PyQt6 | — | GUI 框架（主窗口、分析窗口） |
| pyvista | — | 3D 可视化（探测器几何、粒子径迹） |
| pyvistaqt | — | PyQt6 + pyvista 集成 |
| matplotlib | — | E-θ 二维图和激发能谱绘制 |
| tomli | — | Python < 3.11 的 TOML 解析 |
| numpy | — | 批量模拟中的数组操作和直方图 |
| toml | — | 写入 config 快照（仅批量模式） |

## 13. 常见问题

### Q1: 为什么实验重建激发能可能与设置值有偏差？

实验重建使用探测器测量数据（含能量分辨模糊 1% FWHM for Si, 3% FWHM for CsI，角度重建由网格/击中位置决定），而理论重建使用精确四动量。两者之间的偏差反映了探测器的实际分辨能力。

### Q2: cms_theta 设置为较小值（如 20°）时，产物角度为什么仍可能较大？

`cms_theta` 限制的是质心系角度，实验室系角度由 Lorentz boost 决定。对于逆运动学反应（较重入射粒子打轻靶），质心系前向的粒子在实验室系中可能仍集中在前向小角度范围内，但具体值取决于反应运动学。

### Q3: 粒子未阻停在探测器中怎么办？

在 GUI 模式中，粒子若穿透所有 Si+CsI 探测器仍未阻停，该事件被视为无效（`exited=True`），不会记录到实验分析窗口中。在批量模拟模式中，`process_detector_response` 返回 `None` 跳过该事件。

### Q4: TSpline3 缓存何时需要重新生成？

当以下任一条件变化时，需删除 `assets/range/` 下对应文件以触发重新生成:
- 修改了 `_SPLINE_POINTS` 或能量采样范围
- pycatima 版本更新导致计算结果变化
- 首次运行新粒子-材料组合时会自动生成

### Q5: 为什么配置文件中的 `excitation_C` 键名不变，但不同反应类型对应不同粒子？

为保持配置文件简洁统一，所有反应类型共用 `excitation_C` 键名:
- Type 1: `excitation_C` → 产物 C 激发（等价于母核 A 激发）
- Type 2: `excitation_C` → 产物 C 激发
- Type 3: `excitation_C` → 中间态 C* 激发，`excitation_F` → 产物 F 激发

不存在 `excitation_D` 参数，D 粒子在所有反应类型中始终为基态。