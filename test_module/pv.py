import sys
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QGridLayout,
    QTextEdit, QPushButton, QVBoxLayout, QHBoxLayout,
    QLabel, QSlider
)
from PyQt6.QtCore import Qt

import matplotlib
matplotlib.use("QtAgg")
matplotlib.rcParams["font.family"] = "Times New Roman"
matplotlib.rcParams["mathtext.fontset"] = "stix"
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Scientific Simulation - 2x2 Layout")
        self.resize(1400, 900)

        central = QWidget()
        self.setCentralWidget(central)
        grid = QGridLayout(central)
        grid.setSpacing(6)

        # ---------- Top-Left: Info Panel ----------
        info_panel = QWidget()
        info_layout = QVBoxLayout(info_panel)
        info_layout.setContentsMargins(8, 8, 8, 8)

        title = QLabel("View Status / Operation Log")
        title.setStyleSheet("font-weight: bold; font-size: 14px;")
        info_layout.addWidget(title)

        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setStyleSheet(
            "background-color: #f5f5f5; font-family: Consolas, monospace;"
        )
        info_layout.addWidget(self.info_box)

        self.btn_refresh = QPushButton("Refresh Info")
        self.btn_refresh.clicked.connect(self.update_info)
        info_layout.addWidget(self.btn_refresh)

        # Square wave Fourier series terms
        n_row = QHBoxLayout()
        n_row.addWidget(QLabel("Terms n ="))
        self.n_label = QLabel("5")
        self.n_label.setMinimumWidth(30)
        n_row.addWidget(self.n_label)
        self.n_slider = QSlider(Qt.Orientation.Horizontal)
        self.n_slider.setRange(1, 20)
        self.n_slider.setValue(5)
        self.n_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.n_slider.setTickInterval(1)
        self.n_slider.valueChanged.connect(self.on_draw)
        n_row.addWidget(self.n_slider)
        info_layout.addLayout(n_row)

        # Gaussian sigma parameter
        s_row = QHBoxLayout()
        s_row.addWidget(QLabel("sigma ="))
        self.s_label = QLabel("1.0")
        self.s_label.setMinimumWidth(30)
        s_row.addWidget(self.s_label)
        self.s_slider = QSlider(Qt.Orientation.Horizontal)
        self.s_slider.setRange(1, 30)
        self.s_slider.setValue(10)
        self.s_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.s_slider.setTickInterval(5)
        self.s_slider.valueChanged.connect(self.on_draw)
        s_row.addWidget(self.s_slider)
        info_layout.addLayout(s_row)

        self.btn_draw = QPushButton("Draw")
        self.btn_draw.clicked.connect(self.on_draw)
        self.btn_draw.setStyleSheet("font-weight: bold;")
        info_layout.addWidget(self.btn_draw)

        self.btn_reset = QPushButton("Reset Camera")
        self.btn_reset.clicked.connect(self.reset_cameras)
        info_layout.addWidget(self.btn_reset)

        grid.addWidget(info_panel, 0, 0)

        # ---------- Three Views ----------
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax2d = self.fig.add_subplot(111)
        self.canvas_2d = FigureCanvas(self.fig)  # Top-Right: 2D Fourier
        self.toolbar_2d = NavigationToolbar(self.canvas_2d, self)

        # container: toolbar + canvas
        container_2d = QWidget()
        container_layout = QVBoxLayout(container_2d)
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.setSpacing(0)
        container_layout.addWidget(self.toolbar_2d)
        container_layout.addWidget(self.canvas_2d)

        # Bottom-Right: 3D Grid Surface
        self.plot_gauss3d = QtInteractor(self)

        # Bottom-Left: 1D Histogram
        self.fig_gauss = Figure(figsize=(5, 4), dpi=100)
        self.ax_gauss = self.fig_gauss.add_subplot(111)
        self.canvas_gauss = FigureCanvas(self.fig_gauss)
        self.toolbar_gauss = NavigationToolbar(self.canvas_gauss, self)

        container_gauss = QWidget()
        gauss_layout = QVBoxLayout(container_gauss)
        gauss_layout.setContentsMargins(0, 0, 0, 0)
        gauss_layout.setSpacing(0)
        gauss_layout.addWidget(self.toolbar_gauss)
        gauss_layout.addWidget(self.canvas_gauss)

        grid.addWidget(container_2d, 0, 1)
        grid.addWidget(container_gauss, 1, 0)
        grid.addWidget(self.plot_gauss3d.interactor, 1, 1)

        # row/col stretch
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)

        # draw scenes
        self.draw_surface()
        self.draw_gauss_2d()
        self.draw_gauss_3d()

        # init info
        self.update_info()

    # ---------- View 1: Square Wave Fourier Series (2D) ----------
    def draw_surface(self, n=5):
        x = np.linspace(-np.pi, np.pi, 500)
        fourier = np.zeros_like(x)
        for k in range(1, n + 1, 2):
            fourier += np.sin(k * x) / k
        y = (4.0 / np.pi) * fourier

        self.ax2d.clear()
        self.ax2d.plot(x, y, color="#2c7fb8", linewidth=1.5)
        self.ax2d.set_xlim(-np.pi, np.pi)
        self.ax2d.set_ylim(-1.6, 1.6)
        self.ax2d.set_xlabel("x")
        self.ax2d.set_ylabel("f(x)")
        self.ax2d.set_title(f"Square Wave Fourier Series  n = {n}")
        self.ax2d.grid(True, alpha=0.3)
        self.ax2d.axhline(y=0, color="gray", linewidth=0.5)
        self.ax2d.axvline(x=0, color="gray", linewidth=0.5)
        self.fig.tight_layout()
        self.canvas_2d.draw()

    def on_draw(self):
        n = self.n_slider.value()
        sigma = self.s_slider.value() / 10.0
        self.n_label.setText(str(n))
        self.s_label.setText(f"{sigma:.1f}")
        self.draw_surface(n)
        self.draw_gauss_2d(sigma)
        self.draw_gauss_3d(sigma)
        self.update_info()

    # ---------- View 2: 1D Histogram (TH1D style) ----------
    def draw_gauss_2d(self, sigma=1.0):
        data = np.random.default_rng(42).normal(0, sigma, size=5000)

        self.ax_gauss.clear()
        self.ax_gauss.hist(
            data, bins=50, density=False,
            histtype="step", color="#2c7fb8", linewidth=1.5,
        )
        self.ax_gauss.set_xlabel("Value")
        self.ax_gauss.set_ylabel("Counts")
        self.ax_gauss.set_title(f"1D Histogram (Normal)  sigma = {sigma:.1f}")
        self.ax_gauss.grid(True, alpha=0.3)
        self.fig_gauss.tight_layout()
        self.canvas_gauss.draw()

    # ---------- View 3: Grid Surface (3D) ----------
    def draw_gauss_3d(self, sigma=1.0):
        x, y = np.meshgrid(
            np.linspace(-4, 4, 60),
            np.linspace(-4, 4, 60),
        )
        z = np.sin(x) * np.cos(y)
        grid = pv.StructuredGrid(x, y, z)

        self.plot_gauss3d.clear()
        self.plot_gauss3d.add_mesh(
            grid, cmap="viridis", smooth_shading=True,
            show_scalar_bar=True, scalar_bar_args={"title": "Z Value"},
            show_edges=True, edge_color="#333333",
        )
        self.plot_gauss3d.add_text(
            "Grid Surface  Z = sin(x)*cos(y)",
            position="upper_edge", font_size=10,
        )
        self.plot_gauss3d.enable_anti_aliasing("ssaa")
        self.plot_gauss3d.show_axes()
        self.plot_gauss3d.reset_camera()

    # ---------- Info Panel Update ----------
    def update_info(self):
        sigma = self.s_slider.value()
        lines = []
        lines.append("=== View Overview ===")
        lines.append("[Top-Right] Square Wave Fourier Series (2D)")
        lines.append(f"  Terms: n = {self.n_slider.value()}")
        lines.append("")
        lines.append("[Bottom-Left] 1D Histogram (Normal)")
        lines.append(f"  sigma = {sigma}")
        lines.append("")
        lines.append("[Bottom-Right] Lorenz Attractor (3D)")
        lines.append(f"  sigma = {sigma},  rho = 28,  beta = 8/3")
        lines.append("")
        lines.append("Tip: Drag/scroll in 3D views to rotate/zoom.")
        self.info_box.setPlainText("\n".join(lines))

    def reset_cameras(self):
        self.update_info()


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()