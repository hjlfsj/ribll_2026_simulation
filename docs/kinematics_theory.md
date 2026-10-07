# 运动学理论基础

本文档介绍 RIBLL2026 模拟框架所涉及的相对论运动学基础知识、三种反应类型的生成与抽样方法、以及激发能的运动学重建方法。

---

## 1. 四动量和不变质量基础知识

### 1.1 四动量的定义

在相对论运动学中，粒子的运动状态由**四动量**（four-momentum）统一描述。四动量是一个四维矢量，将能量和三维动量合并：

$$
P = (E, p_x, p_y, p_z) = (E, \vec{p})
$$

采用自然单位制（$c = 1$），四动量的各个分量具有能量量纲（MeV）。在 ROOT 框架中，使用 `TLorentzVector` 类来表示和操作四动量。

### 1.2 四动量的基本性质

**能量-动量关系**（在壳条件，on-shell condition）：

$$
E^2 = p^2 + m^2
$$

其中 $p = |\vec{p}|$ 为动量大小，$m$ 为粒子的静质量（基态质量）。粒子的动能 $E_k$ 与总能量 $E$ 的关系为：

$$
E_k = E - m
$$

**四动量加法**：对于多个粒子组成的系统，系统的总四动量为各粒子四动量之和：

$$
P_{\text{total}} = \sum_i P_i = \left(\sum_i E_i,\; \sum_i \vec{p}_i\right)
$$

### 1.3 不变质量

一个粒子系统（或单个粒子）的**不变质量**（invariant mass）定义为系统总四动量的 Minkowski 模：

$$
M = \sqrt{E_{\text{total}}^2 - p_{\text{total}}^2}
$$

不变质量是 Lorentz 不变量——在所有惯性参考系中数值相同。这是其最重要的性质。

对于单个在壳粒子：

$$
M = \sqrt{E^2 - p^2} = m
$$

即不变质量等于其静质量。对于由多个粒子组成的系统，不变质量反映了系统整体的"有效静质量"。

**激发态粒子**：当粒子处于激发态时，其实际质量为：

$$
m^* = m_{gs} + E_x
$$

其中 $m_{gs}$ 为基态质量，$E_x$ 为激发能（excitation energy）。

### 1.4 Lorentz 变换与 Boost

从一个惯性参考系 $S$ 变换到以速度 $\vec{\beta}$ 相对于 $S$ 运动的参考系 $S'$，四动量的变换规则为 Lorentz boost：

$$
\begin{aligned}
E' &= \gamma (E - \vec{\beta} \cdot \vec{p}) \\
\vec{p}'_{\parallel} &= \gamma (\vec{p}_{\parallel} - \vec{\beta} E) \\
\vec{p}'_{\perp} &= \vec{p}_{\perp}
\end{aligned}
$$

其中 $\gamma = 1 / \sqrt{1 - \beta^2}$。在核反应运动学中，最常用的两个参考系是：

- **实验室系**（lab frame）：靶核静止，入射粒子以束流能量运动
- **质心系**（center-of-mass frame, CMS）：系统总动量 $\vec{P}_{\text{total}} = 0$

质心系的 boost 矢量由系统总四动量决定：

$$
\vec{\beta}_{cm} = \frac{\vec{P}_{\text{total}}}{E_{\text{total}}}
$$

### 1.5 质心系总能量与 Mandelstam $s$

两体碰撞 $A + B$ 的 Mandelstam 变量 $s$ 定义为系统总四动量的平方：

$$
s = (P_A + P_B)^2 = M^2(P_A + P_B)
$$

$\sqrt{s}$ 即为质心系总能量。两体反应的阈值条件为：

$$
\sqrt{s} \geq m_C + m_D
$$

即质心系总能量必须不小于产物粒子的静质量之和。

### 1.6 Q 值与反应能量条件

反应的 **Q 值**定义为初态与末态静质量之差：

$$
Q = \sum_i m_i^{\text{initial}} - \sum_f m_f^{\text{final}}
$$

- $Q > 0$：放热反应（exothermic），无能量阈值
- $Q < 0$：吸热反应（endothermic），需满足束流能量大于阈能

对于两体反应 $A + B \to C + D$，当 $Q < 0$ 时，实验室系阈能为：

$$
E_{\text{threshold}} = -Q \cdot \frac{m_A + m_B}{m_B}
$$

### 1.7 质心系动量

对于两体反应 $A + B \to C + D$，在质心系中，两个出射粒子动量大小相等、方向相反。质心系动量大小由两体相空间公式给出：

$$
p_{cm} = \frac{\sqrt{[s - (m_C + m_D)^2]\,[s - (m_C - m_D)^2]}}{2\sqrt{s}}
$$

质心系中各粒子的能量为：

$$
E_C^{cm} = \sqrt{p_{cm}^2 + m_C^2}, \quad E_D^{cm} = \sqrt{p_{cm}^2 + m_D^2}
$$

---

## 2. 三种反应类型的生成与随机抽样

本项目的运动学模拟涵盖三种核反应类型。所有运动学计算均基于相对论四动量（`ROOT.TLorentzVector`）。

| 反应类型 | 反应式 | 模拟方法 |
|:---:|:---|:---|
| Type 1 | $A^* \to B + C$ | `TGenPhaseSpace` 各向同性衰变 |
| Type 2 | $A + B \to C^{(*)} + D$ | 质心系角度采样 + Lorentz boost |
| Type 3 | $A + B \to C^* + D$, $C^* \to E + F$ | Type 2 + `TGenPhaseSpace` 级联 |

---

### 2.1 Type 1: $A^* \to B + C$（衰变）

**反应描述**：母核 $A$ 处于激发态（质量为 $m_A = m_{A,gs} + E_x$），以动能 $E_k$ 沿 $z$ 轴飞行，在飞行过程中衰变为产物 $B$ 和 $C$。$B$ 和 $C$ 始终处于基态。

**抽样方法**：

使用 ROOT 的 `TGenPhaseSpace` 类进行均匀相空间抽样。该方法在母核 $A$ 的静止系中按相空间分布各向同性地生成衰变产物。

**具体步骤**：

1. **构造母核四动量**：
   $$
   P_A = (E_A, 0, 0, p_A), \quad E_A = E_k + m_A, \quad p_A = \sqrt{E_A^2 - m_A^2}
   $$

2. **相空间生成**：调用 `TGenPhaseSpace.SetDecay(P_A, 2, [m_B, m_C])` 设置衰变，然后 `Generate()` 进行抽样。生成的两个粒子的四动量 $P_B$、$P_C$ 直接位于实验室系中。

3. **角度分布**：`TGenPhaseSpace` 在母核静止系中各向同性，无质心系角度限制。

**有效性检查**：
- 激发能必须满足 $E_x > -Q_{gs}$（基态 Q 值的绝对值），否则衰变在能量上不可行
- 等效于要求 $Q = m_A - m_B - m_C > 0$

---

### 2.2 Type 2: $A + B \to C^{(*)} + D$（两体反应）

**反应描述**：入射粒子 $A$ 以动能 $E_k$ 撞击静止靶核 $B$（位于实验室系原点），生成出射粒子 $C$ 和 $D$。根据粒子种类和激发能可分为三类：

| 类别 | 条件 |
|:---|:---|
| 弹性散射 | $C = A$, $D = B$, $E_x(C) = 0$ |
| 非弹性散射 | $C = A$, $D = B$, $E_x(C) > 0$ |
| 转移反应 | $C \neq A$ 或 $D \neq B$ |

其中仅 $C$ 可设置为激发态（质量为 $m_C = m_{C,gs} + E_x(C)$），$D$ 始终为基态。

**抽样方法**：

采用**质心系角度采样法**。在质心系中按均匀分布随机抽样出射粒子的方向角，然后 boost 回实验室系。此方法比 `TGenPhaseSpace` 更灵活，支持质心系角度限制。

**具体步骤**：

1. **构造实验室系初态四动量**：
   $$
   P_A = (E_A, 0, 0, p_A), \quad P_B = (m_B, 0, 0, 0)
   $$

2. **计算质心系运动学量**：
   $$
   \begin{aligned}
   \sqrt{s} &= M(P_A + P_B) \\
   p_{cm} &= \frac{\sqrt{[s - (m_C + m_D)^2]\,[s - (m_C - m_D)^2]}}{2\sqrt{s}}
   \end{aligned}
   $$

3. **质心系角度均匀抽样**：
   - $\theta_{cm}$ 在 $(0, \theta_{cm}^{\text{max}})$ 内均匀抽样（由 `cms_theta` 参数控制，默认 $180^\circ$ 全覆盖）
   - $\phi_{cm}$ 在 $(0, 2\pi)$ 内均匀抽样

4. **构建质心系四动量**：
   $$
   \begin{aligned}
   P_C^{cm} &= (E_C^{cm},\; p_{cm} \sin\theta_{cm} \cos\phi_{cm},\; p_{cm} \sin\theta_{cm} \sin\phi_{cm},\; p_{cm} \cos\theta_{cm}) \\
   P_D^{cm} &= (E_D^{cm},\; -p_{cm} \sin\theta_{cm} \cos\phi_{cm},\; -p_{cm} \sin\theta_{cm} \sin\phi_{cm},\; -p_{cm} \cos\theta_{cm})
   \end{aligned}
   $$
   其中 $E_i^{cm} = \sqrt{p_{cm}^2 + m_i^2}$。

5. **Boost 回实验室系**：
   $$
   P_C^{lab} = \text{Boost}(P_C^{cm}, \vec{\beta}_{cm}), \quad P_D^{lab} = \text{Boost}(P_D^{cm}, \vec{\beta}_{cm})
   $$
   其中 $\vec{\beta}_{cm} = \vec{P}_{\text{total}} / E_{\text{total}}$ 为质心系 boost 矢量。

**有效性检查**：
- 若 $Q < 0$，需确保 $E_k > E_{\text{threshold}}$
- $\sqrt{s} \geq m_C + m_D$

**角度限制说明**：`cms_theta` 限制的是产物 $C$ 在 $A+B$ **质心系**中的极角范围。在逆运动学反应（较重入射粒子打轻靶）中，质心系前向角经 Lorentz boost 后在实验室系中仍然集中在前向小角度范围，但具体映射关系取决于反应运动学。

---

### 2.3 Type 3: $A + B \to C^* + D$, $C^* \to E + F$（级联反应）

**反应描述**：分为两步——

- **第一步**：$A + B \to C^* + D$，与 Type 2 完全一致，生成激发态 $C^*$ 和基态 $D$。$C^*$ 的质量为 $m_{C^*} = m_{C,gs} + E_x(C)$。
- **第二步**：$C^* \to E + F$，激发态 $C^*$ 衰变为 $E$（始终基态）和 $F$（可设激发态 $E_x(F)$）。

**抽样方法**：

第一步采用质心系角度采样法（同 Type 2），支持 `cms_theta` 角度限制；第二步使用 `TGenPhaseSpace` 在 $C^*$ 静止系中进行各向同性衰变。

**具体步骤**：

1. **第一步**：调用 Type 2 的模拟函数，生成 $P_{C^*}$ 和 $P_D$。

2. **第二步**：以 $P_{C^*}$ 为母核四动量，调用 `TGenPhaseSpace.SetDecay(P_Cstar, 2, [m_E, m_F])`，其中：
   $$
   m_E = m_{E,gs}, \quad m_F = m_{F,gs} + E_x(F)
   $$

3. `TGenPhaseSpace` 自动在 $C^*$ 的静止系中生成 $P_E$ 和 $P_F$，并直接返回实验室系四动量（该类已内置了从母核静止系到实验室系的 boost 变换）。

4. 第二步衰变的 Q 值为：
   $$
   Q_{\text{breakup}} = m_{C^*} - (m_E + m_F) = m_{C^*} - m_{E,gs} - m_{F,gs} - E_x(F)
   $$

**有效性检查**：
- 第一步需满足 Type 2 的所有条件
- 第二步需满足 $Q_{\text{breakup}} \geq 0$

---

## 3. 运动学重建方法

在实验中，探测器测量的是出射粒子的动能和出射角度。运动学重建的目标是从这些观测量出发，反推出反应中激发态粒子的激发能。以下是三种反应类型各自的重建方法。

> **关于 $\gamma$ 光子动量的近似**：激发态核退激时放出的 $\gamma$ 光子动量量级为 $p_\gamma \sim E_x / c \sim \mathcal{O}(1\,\text{MeV}/c)$，远小于重离子动量 $p \sim \mathcal{O}(10^2\,\text{MeV}/c)$。在以下所有重建方法中，忽略 $\gamma$ 光子动量所引入的相对误差 $\ll 1\%$。

---

### 3.1 Type 1: 直接不变质量法

**反应**：$A^* \to B + C$

**已知量**：探测器测量到 $B$ 和 $C$ 的动能 $E_k(B)$、$E_k(C)$，以及出射角度 $\theta(B)$、$\phi(B)$、$\theta(C)$、$\phi(C)$。

**重建步骤**：

1. 利用测量到的动能和角度，重构 $B$ 和 $C$ 的四动量：
   $$
   \begin{aligned}
   P_B &= (E_B,\; \vec{p}_B), \quad E_B = E_k(B) + m_{B,gs} \\
   P_C &= (E_C,\; \vec{p}_C), \quad E_C = E_k(C) + m_{C,gs}
   \end{aligned}
   $$

2. 计算 $B+C$ 系统的不变质量，即母核 $A$ 的重建质量：
   $$
   M_A^{\text{reco}} = M(P_B + P_C) = \sqrt{(E_B + E_C)^2 - |\vec{p}_B + \vec{p}_C|^2}
   $$

3. 母核 $A$ 的激发能：
   $$
   E_x(A) = M_A^{\text{reco}} - m_{A,gs}
   $$

此方法无需束流能量信息，仅依赖对两个衰变产物的完整运动学测量。

---

### 3.2 Type 2: 两种互补方法

**反应**：$A + B \to C^* + D \quad [C^* \to C_{gs} + \gamma]$

**已知量**：探测器测量到 $C$ 和 $D$ 的动能及出射角度。$D$ 处于基态。

#### 方法一：利用束流能量——缺失质量法

利用已知的束流能量和测量到的 $D$ 粒子信息，通过动量/能量守恒反推 $C^*$ 的四动量。

**重建步骤**：

1. 由束流能量构造入射粒子 $A$ 的四动量：
   $$
   P_A = (E_A, 0, 0, p_A), \quad E_A = E_{\text{beam}} + m_{A,gs}
   $$

2. 由测量数据重构 $D$ 的四动量：
   $$
   P_D = (E_D,\; \vec{p}_D), \quad E_D = E_k(D) + m_{D,gs}
   $$

3. 靶核 $B$ 静止于实验室系：
   $$
   P_B = (m_{B,gs}, 0, 0, 0)
   $$

4. 缺失质量法求得 $C^*$ 的四动量：
   $$
   P_{C^*} = P_A + P_B - P_D
   $$

5. $C^*$ 的激发能：
   $$
   E_x(C) = M(P_{C^*}) - m_{C,gs} = \sqrt{E_{C^*}^2 - |\vec{p}_{C^*}|^2} - m_{C,gs}
   $$

**方法一优缺点**：
- **优点**：只需测量 $D$ 粒子，不依赖对 $C$ 的测量
- **缺点**：束流能量存在 $\sim \text{MeV}$ 量级的能散，直接使用 $E_{\text{beam}}$ 会引入系统误差

#### 方法二：不利用束流能量——动量守恒反推 + 不变质量法

利用动量守恒从末态所有粒子的动量反推束流动量，从而重建 $A$ 的能量，绕过束流能散。

**关键思路**：在 $A + B \to C^* + D$ 反应中，忽略 $\gamma$ 光子动量后有 $\vec{p}_A = \vec{p}_C + \vec{p}_D$。由此可以从 $C$ 和 $D$ 的测量动量反推出入射粒子 $A$ 的动量，再用相对论能量-动量关系得到 $E_A$。

**重建步骤**：

1. 由测量数据重构 $C$ 和 $D$ 的四动量：
   $$
   P_C = (E_C,\; \vec{p}_C), \quad P_D = (E_D,\; \vec{p}_D)
   $$

2. 利用横动量守恒反推 $A$ 的动量（忽略 $\gamma$ 光子动量）：
   $$
   \vec{p}_A^{\text{reco}} = \vec{p}_C + \vec{p}_D
   $$

3. 利用在壳条件重建 $A$ 的总能量：
   $$
   E_A^{\text{reco}} = \sqrt{|\vec{p}_A^{\text{reco}}|^2 + m_{A,gs}^2}
   $$

4. 构造反推的 $A$ 四动量：
   $$
   P_A^{\text{reco}} = (E_A^{\text{reco}},\; \vec{p}_A^{\text{reco}})
   $$

5. 缺失质量法求 $C^*$ 四动量：
   $$
   P_{C^*} = P_A^{\text{reco}} + P_B - P_D
   $$

6. $C^*$ 的激发能：
   $$
   E_x(C) = M(P_{C^*}) - m_{C,gs}
   $$

**方法二优缺点**：
- **优点**：全程不依赖 $E_{\text{beam}}$ 的数值，完全规避束流能散引入的系统误差
- **缺点**：需要同时测量 $C$ 和 $D$ 两个粒子；忽略 $\gamma$ 光子动量会引入极小误差（$\ll 1\%$）

> **本项目代码实现默认采用方法二**，即动量守恒反推束流能量的缺失质量法。这使得模拟的激发能重建精度不受束流能散影响。

---

### 3.3 Type 3: 分场景重建 $E_x(C)$ 和 $E_x(F)$

**反应**：$A + B \to C^* + D$, $C^* \to E + F$

**重建目标**：同时得到 $C$ 的激发能 $E_x(C)$ 和 $F$ 的激发能 $E_x(F)$。

**关键物理图像**：

- $C^*$（质量为 $m_{C^*} = m_{C,gs} + E_x(C)$）通过衰变 $C^* \to E + F$ **释放了全部激发能**，衰变产物中 $E$ 始终处于基态。
- 只有 $F$ 可能处于激发态（$m_F = m_{F,gs} + E_x(F)$），随后 $F^* \to F_{gs} + \gamma$ 放出 $\gamma$ 光子退激。因此整个反应链中，$\gamma$ 光子仅来自 $F$ 的退激。
- $C^*$ 的质量由 $E$ 和 $F$ 的静质量完全确定：
  $$
  m_{C^*} = m_{E,gs} + m_F = m_{E,gs} + m_{F,gs} + E_x(F)
  $$

**重建分为两步**：

| 步骤 | 目标 | 方法 | 公式 |
|:---:|:---|:---|:---|
| 第一步 | $E_x(C)$ | 不变质量法 | $E_x(C) = M(P_E + P_F) - m_{C,gs}$ |
| 第二步 | $E_x(F)$ | $Q_0 - Q_1$ 关系 | $E_x(F) = Q_0 - Q_1$ |

其中：
- $Q_0 = (m_{A,gs} + m_{B,gs}) - (m_{E,gs} + m_{F,gs} + m_D)$：末态粒子**全部基态**时的总反应 Q 值（从质量表直接查表得到）
- $Q_1 = (m_{A,gs} + m_{B,gs}) - (m_{C^*} + m_D)$：第一步反应 Q 值（需实验测定）

**$E_x(F) = Q_0 - Q_1$ 的推导**：

将 $m_{C^*} = m_{E,gs} + m_{F,gs} + E_x(F)$ 代入 $Q_1$ 表达式：

$$
\begin{aligned}
Q_0 - Q_1 &= \big[(m_A+m_B) - (m_{E,gs}+m_{F,gs}+m_D)\big] - \big[(m_A+m_B) - (m_{C^*}+m_D)\big] \\
          &= m_{C^*} - m_{E,gs} - m_{F,gs} \\
          &= E_x(F)
\end{aligned}
$$

等价地：

$$
\boxed{E_x(F) = M(P_E + P_F) - m_{E,gs} - m_{F,gs}}
$$

两步的共同输入都是 $M(P_E + P_F)$（即 $C^*$ 的不变质量），差别在于第二步还需要 $Q_1$。$Q_1$ 的获取方式取决于是否测量到了轻粒子 $D$，下面分两种场景讨论。

> **关于 $\gamma$ 光子动量**：唯一需要近似掉 $\gamma$ 动量的地方是计算 $M(P_E + P_F)$——$F$ 退激放出的 $\gamma$ 动量 $p_\gamma \sim E_x(F)/c \sim \mathcal{O}(1\,\text{MeV}/c)$ 远小于 $p_E, p_F$，忽略其动量引入的误差 $\ll 1\%$。

---

#### 场景一：未测量到轻粒子 $D$（仅测量 $E$ 和 $F$）

当 $D$ 出射角度过大或能量过低而未被探测器阵列覆盖时，$E_x(C)$ 可照常得到，但 $Q_1$ 需要借助束流能量。

**重建步骤**：

1. 由 $E$ 和 $F$ 的测量数据构建四动量（$\gamma$ 动量近似为零）：
   $$
   P_E = (E_E,\; \vec{p}_E), \quad P_F = (E_F,\; \vec{p}_F)
   $$

2. **得到 $C$ 的激发能**（不变质量法，无需束流能量）：
   $$
   E_x(C) = M(P_E + P_F) - m_{C,gs}
   $$

3. 由束流能量构造 $A$ 的四动量：
   $$
   P_A = (E_{\text{beam}} + m_{A,gs},\; 0, 0, p_A), \quad P_B = (m_{B,gs}, 0, 0, 0)
   $$

4. 由第一步反应反推 $D$ 的四动量：
   $$
   P_D^{\text{reco}} = P_A + P_B - P_{C^*}, \quad P_{C^*} = P_E + P_F
   $$

5. 由反推出的 $D$ 得到第一步反应 Q 值。由于 $D$ 始终处于基态，其不变质量应等于基态质量：
   $$
   M(P_D^{\text{reco}}) = m_{D,gs}
   $$
   $$
   Q_1 = (m_{A,gs} + m_{B,gs}) - \big(M(P_E + P_F) + M(P_D^{\text{reco}})\big)
   $$

6. **得到 $F$ 的激发能**：
   $$
   E_x(F) = Q_0 - Q_1
   $$
   其中 $Q_0 = (m_{A,gs} + m_{B,gs}) - (m_{E,gs} + m_{F,gs} + m_{D,gs})$ 由质量表给出。

---

#### 场景二：测量到了轻粒子 $D$（测量 $E$、$F$ 和 $D$）

当 $D$ 也在探测器接受度内并被成功测量时，$Q_1$ 可完全不依赖束流能量。

**重建步骤**：

1. 由 $E$ 和 $F$ 的测量数据构建四动量（同场景一）。

2. **得到 $C$ 的激发能**（不变质量法）：
   $$
   E_x(C) = M(P_E + P_F) - m_{C,gs}
   $$

3. 由 $D$ 的测量数据构建 $P_D$。由于 $D$ 处于基态：
   $$
   M(P_D) = m_{D,gs}
   $$

4. 利用动量守恒反推 $A$ 的四动量（不依赖 $E_{\text{beam}}$）：
   $$
   \vec{p}_A^{\text{reco}} = \vec{p}_E + \vec{p}_F + \vec{p}_D
   $$
   $$
   E_A^{\text{reco}} = \sqrt{|\vec{p}_A^{\text{reco}}|^2 + m_{A,gs}^2}
   $$

5. 由第一步反应得到 $Q_1$（无需束流能量）：
   $$
   Q_1 = (m_{A,gs} + m_{B,gs}) - \big(M(P_E + P_F) + M(P_D)\big)
   $$

6. **得到 $F$ 的激发能**：
   $$
   E_x(F) = Q_0 - Q_1
   $$

**场景二优缺点**：
- **优点**：不依赖束流能量，规避能散误差；可利用 $D$ 的测量进行交叉检验
- **缺点**：要求 $D$ 也在探测器接受度内

---

### 3.4 重建方法总结

| 反应类型 | 重建目标 | 主要方法 | 是否需要 $E_{\text{beam}}$ | $\gamma$ 近似 |
|:---:|:---|:---|:---:|:---|
| Type 1 | $E_x(A)$ | 不变质量法 $M(B+C) - m_{A,gs}$ | 否 | 无 $\gamma$ |
| Type 2 | $E_x(C)$ | 缺失质量法（动量守恒反推束流） | 否（方法二） | 忽略 $C^* \to C_{gs}+\gamma$ 的 $\gamma$ 动量 |
| Type 3（缺 $D$） | $E_x(C)$ + $E_x(F)$ | ①不变质量法 $\to E_x(C)$，②$Q_0 - Q_1 \to E_x(F)$（$Q_1$ 需 $E_{\text{beam}}$） | **是**（仅第②步） | 忽略 $F^* \to F_{gs}+\gamma$ 的 $\gamma$ 动量 |
| Type 3（有 $D$） | $E_x(C)$ + $E_x(F)$ | ①不变质量法 $\to E_x(C)$，②$Q_0 - Q_1 \to E_x(F)$（动量守恒得 $Q_1$） | 否 | 忽略 $F^* \to F_{gs}+\gamma$ 的 $\gamma$ 动量 |

> **核心原则总结**：
>
> - **Type 1** 最为直接：测量两个衰变产物，不变质量即得激发能。
> - **Type 2** 的关键技巧是用 $\vec{p}_C + \vec{p}_D$ 反推 $\vec{p}_A$，从而避开束流能散（$\sim \text{MeV}$）引入的系统误差。$C^*$ 退激放出的 $\gamma$ 动量被忽略。
> - **Type 3** 需同时重建两个激发能：第一步用不变质量法得 $E_x(C) = M(P_E+P_F) - m_{C,gs}$；第二步利用 $E_x(F) = Q_0 - Q_1$ 得 $F$ 的激发能，其中 $Q_0$ 由质量表给出、$Q_1$ 由实验测定。$C^*$ 已通过衰变释放全部激发能，$\gamma$ 仅来自 $F^* \to F_{gs} + \gamma$ 且其动量被忽略。
> - 所有涉及 $\gamma$ 退激的场景中，$p_\gamma \sim \mathcal{O}(1\,\text{MeV}/c) \ll p_{\text{particle}} \sim \mathcal{O}(10^2\,\text{MeV}/c)$，忽略引入的误差 $\ll 1\%$。

---

## 附录：ROOT 框架中的关键工具

本项目的运动学计算依赖 ROOT 框架中的以下类和函数：

| ROOT 工具 | 用途 |
|:---|:---|
| `TLorentzVector` | 四动量的存储、加减、Lorentz boost、不变质量计算（`.M()`） |
| `TGenPhaseSpace` | $n$ 体衰变的均匀相空间抽样（Type 1、Type 3 第二步） |
| `TLorentzVector::Boost()` | 将四动量从当前参考系变换到另一参考系 |
| `TLorentzVector::BoostVector()` | 返回质心系相对于当前参考系的 boost 矢量 $\vec{\beta}$ |

**ROOT Python 示例**：

```python
import ROOT

# 构造四动量
p4 = ROOT.TLorentzVector(px, py, pz, E)

# 不变质量
M = p4.M()

# 两四动量相加
p4_sum = p4_1 + p4_2

# Boost 到另一参考系
beta = p4_sum.BoostVector()
p4_1.Boost(beta)
```