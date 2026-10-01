"""
RIBLL2026 GUI模块 (PyQt6 + pyvista + matplotlib)
三个独立窗口：3D可视化、理论E-theta图、探测E-theta图
"""

import sys
import os
import time
import random
import numpy as np
from math import sqrt, atan2, pi

import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QGridLayout, QHBoxLayout, QVBoxLayout,
    QTextEdit, QPushButton, QLabel, QCheckBox, QLineEdit, QSplitter,
    QMessageBox,
)
from PyQt6.QtCore import Qt

import matplotlib
matplotlib.use("QtAgg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["mathtext.fontset"] = "stix"
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

from ..data.nuclear_data import get_particle_mass, parse_particle
from ..kinematics import (
    simulate_type1_decay, simulate_type2_reaction, simulate_type3_sequential,
)
from ..kinematics.reconstruction import reconstruct_excitation, reconstruct_excitation_experimental
from ..detector.geometry import Target, DSSD, SSD, DetectorArray
from ..energy_loss.catima_wrapper import calculate_energy_loss


TRACK_COLORS = ['#E74C3C', '#3498DB', '#27AE60', '#9B59B6', '#E67E22']
# 低饱和度探测器层颜色（用于信息面板）
LAYER_COLORS = {
    'C靶/2':  '#C8A96E',  # 低饱和度金
    'DSSD1':  '#9DC3E6',  # 低饱和度蓝
    'DSSD2':  '#A9D18E',  # 低饱和度绿
    'DSSD3':  '#B4A0D4',  # 低饱和度紫
    'DSSD4':  '#D4A0A0',  # 低饱和度红
    'SSD':    '#A0D4D0',  # 低饱和度青
    'CsI':    '#D4C0A0',  # 低饱和度橙
}
DETECTOR_COLOR = '#7F8C8D'
TARGET_COLOR = '#F1C40F'
GRID_COLOR = '#BDC3C7'


def parse_config_simple(config):
    r = config.get('reaction', {})
    d = config.get('detectors', {})
    s = config.get('simulation', {})
    return {
        'reaction_type': r.get('type', 2),
        'particle_A': r.get('particle_A', '14O'),
        'particle_B': r.get('particle_B', '2H'),
        'particle_C': r.get('particle_C', '10C'),
        'particle_D': r.get('particle_D', '6Li'),
        'particle_E': r.get('particle_E', '6He'),
        'particle_F': r.get('particle_F', '4He'),
        'E_beam': r.get('E_beam', 490.0),
        'excitation_C': r.get('excitation_C', 0.0),
        'excitation_F': r.get('excitation_F', 0.0),
        'cms_theta': r.get('cms_theta', 180.0),
        'n_events': s.get('n_events', 10000),
        'target_thickness': d.get('target_thickness_um', 100.0),
        'ppac_sigma': d.get('ppac_sigma_xy', 0.0),
        'si_resolution': d.get('si_resolution_fwhm', 0.01),
    }


def get_product_names(params):
    rt = params['reaction_type']
    if rt == 1:
        return [params['particle_B'], params['particle_C']]
    elif rt == 2:
        return [params['particle_C'], params['particle_D']]
    elif rt == 3:
        return [params['particle_E'], params['particle_F'], params['particle_D']]
    return []


def get_excited_particle_name(params):
    rt = params['reaction_type']
    if rt == 1:
        return params['particle_A']
    elif rt == 2:
        return params['particle_C']
    elif rt == 3:
        return params['particle_C']
    return '?'


class EThetaCollector:

    def __init__(self, name, n_theta=90, n_energy=100, theta_max=180, energy_max=600):
        self.name = name
        self.n_theta = n_theta
        self.n_energy = n_energy
        self.theta_max = theta_max
        self.energy_max = energy_max
        self.hist = np.zeros((n_energy, n_theta), dtype=np.float64)

    def fill(self, theta, energy):
        it = int(theta / self.theta_max * self.n_theta)
        ie = int(energy / self.energy_max * self.n_energy)
        if 0 <= it < self.n_theta and 0 <= ie < self.n_energy:
            self.hist[ie, it] += 1

    def reset(self):
        self.hist.fill(0)


class ExcitationCollector:

    def __init__(self, name, n_bins=100, x_min=0, x_max=100):
        self.name = name
        self.n_bins = n_bins
        self.x_min = x_min
        self.x_max = x_max
        self.values = []

    def fill(self, value):
        self.values.append(value)

    def reset(self):
        self.values.clear()


class AnalysisWindow(QMainWindow):
    """独立的分析画布窗口 (Canvas 2 / Canvas 3)"""

    def __init__(self, title, params, product_names, is_detected=False, parent=None):
        super().__init__(parent)
        self.params = params
        self.product_names = product_names
        self.is_detected = is_detected
        self.n_products = len(product_names)

        self.setWindowTitle(title)
        self.resize(1200, 800)

        self.collectors = []
        for pn in product_names:
            self.collectors.append(EThetaCollector(
                pn, energy_max=params['E_beam'] * 1.2))

        self.ex_collector = ExcitationCollector(
            get_excited_particle_name(params))

        self._init_ui()

    def _init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        # ------ rebin control row ------
        rebin_row = QHBoxLayout()
        rebin_row.addWidget(QLabel("E_x rebin:"))
        rebin_row.addWidget(QLabel("bins"))
        self.rebin_bins = QLineEdit("100")
        self.rebin_bins.setMaximumWidth(50)
        rebin_row.addWidget(self.rebin_bins)
        rebin_row.addWidget(QLabel("min"))
        self.rebin_min = QLineEdit("0")
        self.rebin_min.setMaximumWidth(50)
        rebin_row.addWidget(self.rebin_min)
        rebin_row.addWidget(QLabel("max"))
        self.rebin_max = QLineEdit("100")
        self.rebin_max.setMaximumWidth(50)
        rebin_row.addWidget(self.rebin_max)

        self.btn_rebin = QPushButton("Rebin E_x")
        self.btn_rebin.clicked.connect(self._on_rebin)
        rebin_row.addWidget(self.btn_rebin)
        rebin_row.addStretch()
        main_layout.addLayout(rebin_row)

        # ------ 2x2 matplotlib canvas ------
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)

    def fill_event(self, kin_result, detected_data=None):
        if not self.is_detected:
            particles = self._get_particle_list(kin_result)
            for i, p in enumerate(particles):
                if i < len(self.collectors):
                    self.collectors[i].fill(p['theta'], p['Ek'])
            E_x = reconstruct_excitation(kin_result, self.params)
            self.ex_collector.fill(E_x)
        else:
            if detected_data is not None:
                for i in range(min(len(self.collectors),
                                   detected_data.get('n_particles', 0))):
                    tk = 'det_p{}_theta'.format(i)
                    ek = 'det_p{}_Eexp'.format(i)
                    if tk in detected_data and ek in detected_data:
                        self.collectors[i].fill(detected_data[tk], detected_data[ek])
                if 'Ex_exp' in detected_data:
                    self.ex_collector.fill(detected_data['Ex_exp'])

    def _get_particle_list(self, kin_result):
        rt = self.params['reaction_type']
        if rt == 1:
            return [kin_result['B'], kin_result['C']]
        elif rt == 2:
            return [kin_result['C'], kin_result['D']]
        elif rt == 3:
            return [kin_result['E'], kin_result['F'], kin_result['D']]
        return []

    def reset(self):
        for c in self.collectors:
            c.reset()
        self.ex_collector.reset()

    def draw(self):
        self.fig.clear()

        for i, c in enumerate(self.collectors):
            ax = self.fig.add_subplot(2, 2, i + 1)
            if c.hist.sum() > 0:
                h = c.hist
                h_disp = np.where(h > 0, np.log10(h), -1) if h.max() > 1 else h
                im = ax.pcolormesh(
                    np.linspace(0, c.theta_max, c.n_theta + 1),
                    np.linspace(0, c.energy_max, c.n_energy + 1),
                    h_disp, cmap='viridis', shading='auto',
                )
                self.fig.colorbar(
                    im, ax=ax,
                    label='log10(counts)' if h.max() > 1 else 'counts')
            else:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                        ha='center', va='center', color='gray')
            ax.set_xlabel(r'$\theta_{%s}$ [deg]' % c.name)
            if self.is_detected:
                ax.set_ylabel(r'$E_{%s}^{exp}$ [MeV]' % c.name)
                ax.set_title('{} 实验 (E^{{exp}} vs {})'.format(c.name, r'$\theta$'))
            else:
                ax.set_ylabel(r'$E_{%s}$ [MeV]' % c.name)
                ax.set_title('{} (E vs {})'.format(c.name, r'$\theta$'))

        ax_ex = self.fig.add_subplot(2, 2, 4)
        if len(self.ex_collector.values) > 0:
            try:
                nb = int(self.rebin_bins.text())
                xmin = float(self.rebin_min.text())
                xmax = float(self.rebin_max.text())
            except ValueError:
                nb = self.ex_collector.n_bins
                xmin = self.ex_collector.x_min
                xmax = self.ex_collector.x_max
            ax_ex.hist(
                self.ex_collector.values, bins=nb, range=(xmin, xmax),
                histtype='step', color='#E74C3C', linewidth=1.5,
            )
            ax_ex.set_xlabel(r'$E_x$ [MeV]')
            ax_ex.set_ylabel('Counts')
            prefix = '实验 ' if self.is_detected else '理论 '
            ax_ex.set_title(prefix + r'$E_x(%s)$' % get_excited_particle_name(self.params))
        else:
            ax_ex.text(0.5, 0.5, 'No data', transform=ax_ex.transAxes,
                       ha='center', va='center', color='gray')

        self.fig.tight_layout()
        self.canvas.draw()

    def _on_rebin(self):
        self.draw()


class MainWindow(QMainWindow):
    """主窗口 (Canvas 1): 左侧事件详情 + 右侧3D可视化"""

    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.config = config
        self.params = parse_config_simple(config)
        self.n_simulated = 0

        self.product_names = get_product_names(self.params)

        self._detector_geometry_data = None

        self._init_ui()
        self._draw_3d_geometry()

    def _init_ui(self):
        self.setWindowTitle("RIBLL2026 - 3D Visualization & Event Info")
        self.resize(1400, 800)

        central = QWidget()
        self.setCentralWidget(central)
        vbox = QVBoxLayout(central)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(0)

        # ====== Top bar: Quit button (right-aligned) ======
        top_bar = QHBoxLayout()
        top_bar.addStretch()
        self.btn_quit = QPushButton("Quit")
        self.btn_quit.clicked.connect(QApplication.instance().quit)
        top_bar.addWidget(self.btn_quit)
        vbox.addLayout(top_bar)

        # ====== Main content: Left panel + 3D view ======
        hbox = QHBoxLayout()
        hbox.setSpacing(6)
        vbox.addLayout(hbox)

        # ====== Left panel: Controls + Event Info ======
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(8, 8, 8, 8)

        title = QLabel("RIBLL2026 Controls & Event Info")
        title.setStyleSheet("font-weight: bold; font-size: 14px;")
        left_layout.addWidget(title)

        row = QHBoxLayout()
        self.check_batch = QCheckBox("Batch Mode")
        self.check_batch.setChecked(False)
        row.addWidget(self.check_batch)

        row.addWidget(QLabel("N:"))
        self.edit_n_events = QLineEdit(str(self.params['n_events']))
        self.edit_n_events.setMaximumWidth(70)
        row.addWidget(self.edit_n_events)

        self.btn_sim = QPushButton("Sim")
        self.btn_sim.setStyleSheet("font-weight: bold;")
        row.addWidget(self.btn_sim)

        self.btn_save = QPushButton("Save .root")
        row.addWidget(self.btn_save)
        left_layout.addLayout(row)

        self.label_status = QLabel("Ready")
        self.label_status.setStyleSheet("color: #2980B9; font-weight: bold;")
        left_layout.addWidget(self.label_status)

        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setStyleSheet(
            "background-color: #F8F9FA; font-family: Consolas, monospace; font-size: 12px;")
        left_layout.addWidget(self.info_box)

        hbox.addWidget(left, 1)

        # ====== Right panel: 3D pyvista view ======
        self.plot_3d = QtInteractor(self)
        hbox.addWidget(self.plot_3d.interactor, 2)

    # ===================== 3D Geometry =====================

    def _build_detector_geometry_data(self):
        polys = []
        colors_list = []

        target = Target(radius=15.0, thickness=self.params['target_thickness'])
        half_t = target.thickness_mm / 2.0
        r = target.radius
        h = max(half_t * 2, 0.1)

        cyl = pv.Cylinder(
            center=(0, 0, 0), direction=(0, 0, 1),
            radius=r, height=h, resolution=32,
        )
        polys.append(cyl)
        colors_list.append(TARGET_COLOR)

        array = DetectorArray(
            target=Target(radius=15.0, thickness=self.params['target_thickness']))
        array.build_default()

        for det in array.detectors:
            hw = det.half_w
            hh = det.half_h
            z = det.z
            det_color = LAYER_COLORS.get(det.name, DETECTOR_COLOR)

            plane = pv.Plane(
                center=(0, 0, z), direction=(0, 0, 1),
                i_size=hw * 2, j_size=hh * 2,
                i_resolution=1, j_resolution=1,
            )
            polys.append(plane)
            colors_list.append(det_color)

            if isinstance(det, DSSD):
                gs = det.grid_size
                gn = det.grid_n
                for k in range(1, gn):
                    x = -hw + k * gs
                    polys.append(pv.Line(pointa=(x, -hh, z), pointb=(x, hh, z)))
                    colors_list.append(GRID_COLOR)
                    y = -hh + k * gs
                    polys.append(pv.Line(pointa=(-hw, y, z), pointb=(hw, y, z)))
                    colors_list.append(GRID_COLOR)

        self._detector_geometry_data = (polys, colors_list)

    def _draw_3d_geometry(self):
        self.plot_3d.clear()
        if self._detector_geometry_data is None:
            self._build_detector_geometry_data()
        polys, colors_list = self._detector_geometry_data
        for i, mesh in enumerate(polys):
            self.plot_3d.add_mesh(
                mesh, color=colors_list[i],
                show_edges=(i == 0), edge_color='#555555', line_width=1,
            )
        self.plot_3d.add_text(
            "RIBLL2026 - Detector Setup", position="upper_edge",
            font_size=10, color='black',
        )
        self.plot_3d.show_axes()
        self.plot_3d.view_isometric()
        self.plot_3d.reset_camera()

    def _draw_3d_event(self, kin_result, detected_data=None):
        self.plot_3d.clear()
        if self._detector_geometry_data is None:
            self._build_detector_geometry_data()
        polys, colors_list = self._detector_geometry_data
        for i, mesh in enumerate(polys):
            self.plot_3d.add_mesh(
                mesh, color=colors_list[i],
                show_edges=(i == 0), edge_color='#555555',
                line_width=1, opacity=0.4,
            )

        particles = self._get_particle_list(kin_result)

        # ====== Z-axis (gray, reference line) ======
        array = DetectorArray(
            target=Target(radius=15.0, thickness=self.params['target_thickness']))
        array.build_default()
        z_axis = pv.Line(pointa=(0, 0, -10.0), pointb=(0, 0, 200.0))
        self.plot_3d.add_mesh(
            z_axis, color='#999999', line_width=3, opacity=0.6,
        )

        for i, p in enumerate(particles):
            p4 = p['p4']
            px, py, pz = p4.Px(), p4.Py(), p4.Pz()
            name = p['name']
            color = TRACK_COLORS[i % len(TRACK_COLORS)]

            # determine track endpoint from detected data
            ex, ey, ez = 0.0, 0.0, 0.0
            if detected_data is not None and 'particles_info' in detected_data:
                for pi in detected_data['particles_info']:
                    if pi['name'] == name:
                        ox, oy, oz = pi['track_origin']
                        ex, ey, ez = pi['track_end']
                        break
                else:
                    # fallback: use momentum-scaled
                    scale = 80.0 / max(p4.P(), 1.0)
                    ox, oy, oz = 0.0, 0.0, 0.0
                    ex, ey, ez = px * scale, py * scale, pz * scale
            else:
                ox, oy, oz = 0.0, 0.0, 0.0
                scale = 80.0 / max(p4.P(), 1.0)
                ex, ey, ez = px * scale, py * scale, pz * scale

            line = pv.Line(pointa=(ox, oy, oz), pointb=(ex, ey, ez))
            self.plot_3d.add_mesh(line, color=color, line_width=3, label=name)

            sphere = pv.Sphere(
                center=(ex, ey, ez), radius=2.0,
            )
            self.plot_3d.add_mesh(sphere, color=color)

            # label at z=50 mm on track
            t50 = 50.0 / pz if pz > 1e-30 else 0.0
            label_pos = (ox + px * t50, oy + py * t50, 50.0)
            self.plot_3d.add_point_labels(
                [label_pos], [name],
                font_size=10, point_size=1, point_color=color,
                shape_opacity=0.0,
            )

        self.plot_3d.add_text(
            "Event #{}".format(self.n_simulated),
            position="upper_edge", font_size=10, color='black',
        )
        self.plot_3d.show_axes()
        self.plot_3d.view_isometric()

    # ===================== Info Panel =====================

    def _update_info(self, kin_result, detected_data=None):
        if kin_result is None:
            self.info_box.setHtml("<pre>No valid event.</pre>")
            return
        lines = []
        rt = self.params['reaction_type']
        cms = self.params.get('cms_theta', 180.0)

        type_desc = {
            1: "A → B+C  (衰变, A 激发)",
            2: "A+B → C+D  (转移反应, C 激发)",
            3: "A+B → C*+D → E+F+D  (级联衰变, C+F 激发)",
        }
        lines.append("=== Event #{} ===".format(self.n_simulated))
        lines.append("Type: {}  —  {}".format(rt, type_desc.get(rt, "")))
        lines.append("cms_theta: {:.0f}°  (C in CMS: 0–{:.0f}°)".format(cms, cms))
        lines.append("Beam: {} @ {:.1f} MeV".format(
            self.params['particle_A'], self.params['E_beam']))
        lines.append("Target: {}".format(self.params['particle_B']))
        q_val = kin_result.get('Q_value', 0)
        if rt == 3:
            q_val = kin_result.get('step1', {}).get('Q_value', 0)
        lines.append("Q-value: {:.3f} MeV".format(q_val))
        E_x = reconstruct_excitation(kin_result, self.params)
        lines.append("E_x({}): {:.3f} MeV".format(
            get_excited_particle_name(self.params), E_x))

        rp = detected_data.get('reaction_point', None) if detected_data else None
        if rp:
            lines.append("Reaction point: x={:.3f} mm, y={:.3f} mm, z={:.1f} um".format(
                rp[0], rp[1], rp[2]))
        lines.append("")

        particles = self._get_particle_list(kin_result)
        for p in particles:
            lines.append("{}: Ek={:.2f} MeV, theta={:.2f}deg, phi={:.2f}deg".format(
                p['name'], p['Ek'], p['theta'], p['phi']))

        cms_theta_c = kin_result.get('cms_theta_C', None)
        if cms_theta_c is not None:
            lines.append("C CMS angle: {:.2f}°".format(cms_theta_c))

        if detected_data is None or 'particles_info' not in detected_data:
            self.info_box.setHtml("<pre>{}</pre>".format("\n".join(lines)))
            return

        lines.append("")
        lines.append("=" * 56)

        for pi in detected_data['particles_info']:
            lines.append("")
            lines.append("--- {name} (Z={Z}, A={A})  Ek0={Ek0:.2f} MeV ---".format(
                name=pi['name'], Z=pi['Z'], A=pi['A'], Ek0=pi['Ek_initial']))
            fmt_hdr = "  {:<14s} {:>10s} {:>10s} {:>10s} {:>6s}"
            fmt_row = "  {:<14s} {:>10.3f} {:>10s} {:>10.3f} {:>6s}"
            lines.append(fmt_hdr.format(
                "Layer", "ThE-loss", "DetE-loss", "Remain", "状态"))
            lines.append("  " + "-" * 52)

            for ly in pi['layers']:
                det_str = "{:>10.3f}".format(ly['E_detected']) if ly['E_detected'] is not None else "        --"
                status = "阻停" if ly['stopped'] else "穿过"
                color = LAYER_COLORS.get(ly['name'], '#AAAAAA')
                ly_line = fmt_row.format(
                    ly['name'], ly['E_theory'], det_str,
                    ly['remaining'], status)
                lines.append('<span style="color:{};">{}</span>'.format(color, ly_line))

            total_det = pi['E_total_detected']
            target_loss = pi['layers'][0]['E_theory'] if pi['layers'] else 0.0
            E_exp = total_det + target_loss
            lines.append("  " + "-" * 52)
            lines.append("  {:>14s} {:>10.3f} {:>10.3f}  测量theta={:.2f}°".format(
                "Total(Si+CsI)", pi['E_total_theory'], total_det, pi['theta_reco']))
            lines.append("  {:>14s} {:>10.3f}  (靶校正后粒子能量)".format(
                "E_exp", E_exp))
            if pi['exited']:
                lines.append("  !! 粒子离开探测器阵列, 未被探测 !!")
            else:
                lines.append("  阻停在: {}".format(pi['stopped_in']))

        lines.append("")
        lines.append("=" * 56)

        exc_particle = get_excited_particle_name(self.params)
        E_x_set = self.params.get('excitation_C', 0.0)
        E_x_theory = E_x
        E_x_exp = detected_data.get('Ex_exp', 0.0)

        lines.append("E_x({}) 设置值:   {:.3f} MeV".format(exc_particle, E_x_set))
        lines.append("E_x({}) 理论重建: {:.3f} MeV".format(exc_particle, E_x_theory))
        lines.append("E_x({}) 实验重建: {:.3f} MeV".format(exc_particle, E_x_exp))

        self.info_box.setHtml("<pre>{}</pre>".format("\n".join(lines)))

    # ===================== Event Generation =====================

    def _get_particle_list(self, kin_result):
        rt = self.params['reaction_type']
        if rt == 1:
            return [kin_result['B'], kin_result['C']]
        elif rt == 2:
            return [kin_result['C'], kin_result['D']]
        elif rt == 3:
            return [kin_result['E'], kin_result['F'], kin_result['D']]
        return []

    def _generate_event(self):
        p = self.params
        try:
            if p['reaction_type'] == 1:
                return simulate_type1_decay(
                    p['particle_A'], p['E_beam'],
                    p['particle_B'], p['particle_C'],
                    excitation_C=p['excitation_C'],
                )
            elif p['reaction_type'] == 2:
                return simulate_type2_reaction(
                    p['particle_A'], p['E_beam'], p['particle_B'],
                    p['particle_C'], p['particle_D'],
                    excitation_C=p['excitation_C'],
                    cms_theta=p['cms_theta'],
                )
            elif p['reaction_type'] == 3:
                return simulate_type3_sequential(
                    p['particle_A'], p['E_beam'], p['particle_B'],
                    p['particle_C'], p['particle_D'],
                    p['particle_E'], p['particle_F'],
                    excitation_C=p['excitation_C'],
                    excitation_F=p['excitation_F'],
                    cms_theta=p['cms_theta'],
                )
        except Exception as e:
            print("Event generation error: {}".format(e))
        return None

    def _process_detector(self, kin_result):
        p = self.params
        particles = self._get_particle_list(kin_result)
        if not particles:
            return None

        target_thickness = p['target_thickness']
        si_res = p['si_resolution']
        csi_res = 0.03
        ppac_sigma = p['ppac_sigma']

        array = DetectorArray(
            target=Target(radius=15.0, thickness=target_thickness))
        array.build_default()

        pt = array.target.sample_reaction_point(ppac_sigma, ppac_sigma)
        ox, oy, oz = pt['true']
        ox_obs, oy_obs, z_obs = pt['observed']
        half_target_um = target_thickness / 2.0

        detected = {'n_particles': len(particles)}
        particles_info = []
        all_detected = True

        for i_p, p_info in enumerate(particles):
            p4 = p_info['p4']
            px, py, pz = p4.Px(), p4.Py(), p4.Pz()
            Ek = p_info['Ek']
            name = p_info['name']

            if pz <= 0:
                all_detected = False
                break

            Z, A = parse_particle(name)

            Ek_after_target = calculate_energy_loss(
                Z, A, Ek, 'C', half_target_um)
            target_loss = Ek - Ek_after_target

            layers = [{
                'name': 'C靶/2',
                'material': 'C',
                'thickness_um': half_target_um,
                'E_theory': target_loss,
                'E_detected': None,
                'remaining': Ek_after_target,
                'stopped': False,
            }]

            remaining = Ek_after_target
            stopped = False
            stopped_in = None
            theta_reco = p_info['theta']
            last_hx, last_hy, last_det_z = ox, oy, 0.0  # track endpoint tracking

            for det in array.detectors:
                if not hasattr(det, 'material'):
                    continue

                can, hx, hy = det.can_hit(ox, oy, oz, px, py, pz)
                if not can:
                    continue

                if det.material == 'Si':
                    last_hx, last_hy, last_det_z = hx, hy, det.z
                    E_after = calculate_energy_loss(
                        Z, A, remaining, 'Si', det.thickness_um)
                    deposit = remaining - E_after
                    sigma = si_res * deposit / 2.355
                    deposit_smeared = deposit + random.gauss(0, max(sigma, 0.001))

                    layer_stopped = (E_after <= 0.0)

                    if layer_stopped and isinstance(det, DSSD):
                        cx, cy, gi, gj = det.get_grid_center(hx, hy)
                        dx_grid = cx - ox_obs
                        dy_grid = cy - oy_obs
                        dz_det = det.z
                        theta_reco = atan2(
                            sqrt(dx_grid**2 + dy_grid**2), dz_det) * 180.0 / pi
                    elif layer_stopped:
                        dx = hx - ox_obs
                        dy = hy - oy_obs
                        dz_det = det.z
                        theta_reco = atan2(
                            sqrt(dx**2 + dy**2), dz_det) * 180.0 / pi

                    layers.append({
                        'name': det.name,
                        'material': 'Si',
                        'thickness_um': det.thickness_um,
                        'E_theory': deposit,
                        'E_detected': deposit_smeared,
                        'remaining': max(0.0, E_after),
                        'stopped': layer_stopped,
                    })

                    if layer_stopped:
                        stopped = True
                        stopped_in = det.name
                        remaining = 0.0
                        break
                    remaining = E_after

                elif det.material == 'CsI':
                    last_hx, last_hy, last_det_z = hx, hy, det.z
                    deposit = remaining
                    sigma = csi_res * deposit / 2.355
                    deposit_smeared = deposit + random.gauss(0, max(sigma, 0.001))

                    dx = hx - ox_obs
                    dy = hy - oy_obs
                    dz_det = det.z
                    theta_reco = atan2(
                        sqrt(dx**2 + dy**2), dz_det) * 180.0 / pi

                    layers.append({
                        'name': det.name,
                        'material': 'CsI',
                        'thickness_um': det.thickness_um,
                        'E_theory': deposit,
                        'E_detected': deposit_smeared,
                        'remaining': 0.0,
                        'stopped': True,
                    })
                    stopped = True
                    stopped_in = det.name
                    remaining = 0.0
                    break

            if not stopped:
                all_detected = False
                last_det_z = max(d.z for d in array.detectors if hasattr(d, 'z'))
                far_z = last_det_z + 150.0
                t_far = far_z / pz if pz > 1e-30 else 0.0
                last_hx = ox + px * t_far
                last_hy = oy + py * t_far
                last_det_z = far_z

            total_theory = sum(l['E_theory'] for l in layers if l['name'] != 'C靶/2')
            total_detected = sum(
                l['E_detected'] for l in layers
                if l['name'] != 'C靶/2' and l['E_detected'] is not None
            )

            particles_info.append({
                'name': name,
                'Z': Z,
                'A': A,
                'Ek_initial': Ek,
                'theta_reco': theta_reco,
                'E_total_theory': total_theory,
                'E_total_detected': total_detected,
                'layers': layers,
                'stopped_in': stopped_in,
                'exited': not stopped,
                'track_origin': (ox, oy, 0.0),
                'track_end': (last_hx, last_hy, last_det_z),
            })

            detected['det_p{}_theta'.format(i_p)] = theta_reco
            detected['det_p{}_Ek'.format(i_p)] = max(0.0, total_detected)
            detected['det_p{}_Eexp'.format(i_p)] = max(0.0, total_detected + target_loss)

        detected['particles_info'] = particles_info
        detected['all_detected'] = all_detected
        detected['reaction_point'] = (ox, oy, oz * 1000.0)  # mm, mm, um
        detected['reconstructed_E_x'] = reconstruct_excitation(kin_result, p)
        detected['Ex_exp'] = reconstruct_excitation_experimental(detected, kin_result, p)
        return detected

    # ===================== Sim / Save =====================

    def setup_callbacks(self, on_sim_callback, on_save_callback):
        self.btn_sim.clicked.connect(on_sim_callback)
        self.btn_save.clicked.connect(on_save_callback)

    def on_save(self, theory_win, detected_win):
        from datetime import datetime
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = self.config.get('output', {}).get('data_dir', './output')
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, 'sim_{}.root'.format(timestamp))

        import ROOT
        f = ROOT.TFile(path, 'RECREATE')

        for win, prefix in [(theory_win, 'theory'), (detected_win, 'detected')]:
            for i, c in enumerate(win.collectors):
                h = ROOT.TH2D(
                    '{}_{}_p{}'.format(prefix, 'h2', i),
                    '{} (E vs theta) [{}];theta [deg];E [MeV]'.format(
                        c.name, 'Theory' if prefix == 'theory' else 'Detected'),
                    c.n_theta, 0, c.theta_max,
                    c.n_energy, 0, c.energy_max,
                )
                for ie in range(c.n_energy):
                    for it in range(c.n_theta):
                        if c.hist[ie, it] > 0:
                            h.SetBinContent(it + 1, ie + 1, c.hist[ie, it])
                h.Write()

            if win.ex_collector.values:
                h_ex = ROOT.TH1D(
                    '{}_h_ex'.format(prefix),
                    '{} Reconstructed E_x;E_x [MeV];Counts'.format(
                        'Theory' if prefix == 'theory' else 'Detected'),
                    win.ex_collector.n_bins,
                    win.ex_collector.x_min,
                    win.ex_collector.x_max,
                )
                for v in win.ex_collector.values:
                    h_ex.Fill(v)
                h_ex.Write()

        f.Close()
        self.label_status.setText("Saved: {}".format(path))
        QApplication.processEvents()


class RIBLLApp:
    """应用控制器：创建三个窗口并管理模拟逻辑"""

    def __init__(self, config):
        self.config = config
        self.params = parse_config_simple(config)

        self.main_win = MainWindow(config)
        self.main_win.setup_callbacks(self.on_sim, self._on_save)

        self.theory_win = AnalysisWindow(
            "RIBLL2026 - Theory Analysis (E-theta & E_x)",
            self.params,
            get_product_names(self.params),
            is_detected=False,
        )
        self.theory_win.setWindowTitle("RIBLL2026 - Theory Analysis")

        self.detected_win = AnalysisWindow(
            "RIBLL2026 - Detected Analysis (E-theta & E_x)",
            self.params,
            get_product_names(self.params),
            is_detected=True,
        )
        self.detected_win.setWindowTitle("RIBLL2026 - Detected Analysis")

        self._position_windows()
        self.main_win.show()
        self.theory_win.show()
        self.detected_win.show()

    def _position_windows(self):
        screen = QApplication.primaryScreen()
        if screen is None:
            return
        geo = screen.availableGeometry()
        sw = geo.width()
        sh = geo.height()

        self.main_win.resize(int(sw * 0.50), int(sh * 0.60))
        self.main_win.move(int(sw * 0.01), int(sh * 0.02))

        self.theory_win.resize(int(sw * 0.48), int(sh * 0.52))
        self.theory_win.move(int(sw * 0.01), int(sh * 0.62))

        self.detected_win.resize(int(sw * 0.48), int(sh * 0.52))
        self.detected_win.move(int(sw * 0.51), int(sh * 0.62))

    def on_sim(self):
        is_batch = self.main_win.check_batch.isChecked()
        if is_batch:
            self._run_batch()
        else:
            self._run_single()

    def _run_single(self):
        mw = self.main_win
        mw.n_simulated += 1
        mw.label_status.setText("Simulating event {}...".format(mw.n_simulated))
        QApplication.processEvents()

        kin_result = mw._generate_event()
        if kin_result is None or not kin_result.get('valid'):
            mw.label_status.setText("Event {}: INVALID".format(mw.n_simulated))
            mw._update_info(None)
            return

        detected_data = mw._process_detector(kin_result)
        mw._draw_3d_event(kin_result, detected_data)
        mw._update_info(kin_result, detected_data)

        self.theory_win.fill_event(kin_result)
        self.theory_win.draw()

        if detected_data is not None and detected_data.get('all_detected', False):
            self.detected_win.fill_event(kin_result, detected_data)
            self.detected_win.draw()
            mw.label_status.setText("Event {}: OK (detected)".format(mw.n_simulated))
        else:
            mw.label_status.setText("Event {}: OK (partial/not detected)".format(mw.n_simulated))

        QApplication.processEvents()

    def _run_batch(self):
        mw = self.main_win
        try:
            n = int(mw.edit_n_events.text())
        except ValueError:
            n = self.params['n_events']
        self.params['n_events'] = n
        mw.label_status.setText("Batch: 0/{}...".format(n))
        QApplication.processEvents()

        self.theory_win.reset()
        self.detected_win.reset()

        n_detected = 0
        t0 = time.time()
        last_event = None
        last_detected_data = None

        for i in range(n):
            if i % max(1, n // 10) == 0:
                mw.label_status.setText(
                    "Batch: {}/{} ({} detected)".format(i, n, n_detected))
                QApplication.processEvents()

            kin_result = mw._generate_event()
            if kin_result is None or not kin_result.get('valid'):
                continue

            self.theory_win.fill_event(kin_result)

            detected_data = mw._process_detector(kin_result)
            if detected_data is not None and detected_data.get('all_detected', False):
                self.detected_win.fill_event(kin_result, detected_data)
                n_detected += 1
                last_event = kin_result
                last_detected_data = detected_data

        elapsed = time.time() - t0

        self.theory_win.draw()
        self.detected_win.draw()

        if last_event is not None:
            mw._draw_3d_event(last_event, last_detected_data)
            mw._update_info(last_event, last_detected_data)

        rate = n / elapsed if elapsed > 0 else 0
        mw.label_status.setText(
            "Done: {} events in {:.1f}s ({:.0f} evt/s), {} detected ({:.1f}%)".format(
                n, elapsed, rate, n_detected,
                100.0 * n_detected / max(n, 1)))
        QApplication.processEvents()

    def _on_save(self):
        self.main_win.on_save(self.theory_win, self.detected_win)


def launch_gui(config):
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    controller = RIBLLApp(config)
    sys.exit(app.exec())