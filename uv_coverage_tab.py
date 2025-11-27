# uv_coverage_tab.py

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QComboBox, QPushButton
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

import GUI_master_code as ms


class UVCoverageTab(QWidget):
    """
    Tab for displaying the U-V plane coverage of the interferometer.
    """
    def __init__(self, parent=None):
        super().__init__(parent)

        self._have_plot = False

        # Standard baseline colors matching other tabs
        self.baseline_color_list = [
            "#5ac4ee", "#63e3bf", "#b95cf4",
            "#8cf55a", "#f559b7", "#f7b54c"
        ]

        self._create_widgets()
        self._create_layouts()

    def _create_widgets(self):
        # --- Data Selection Group -----------------------------------------
        self.data_group = QGroupBox("Data Selection")

        self.epoch_label = QLabel("Epoch:")
        self.epoch_combo = QComboBox()
        self.epoch_combo.addItems([
            "2018_red", "2019_red", "2020_red",
            "2021_red", "2022a_red", "2022b_red", "2023_red"
        ])

        self.file_label = QLabel("File:")
        self.file_combo = QComboBox()

        self.update_button = QPushButton("Update Plot")
        self.update_button.clicked.connect(self.update_plot)

        # --- Matplotlib Canvas --------------------------------------------
        self.figure = Figure(figsize=(5, 5)) # Square aspect ratio works best for UV plots
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)

        # --- Connections --------------------------------------------------
        self.epoch_combo.currentTextChanged.connect(self.update_files)
        
        # Initialize
        self.update_files()

    def _create_layouts(self):
        # --- Left Side Layout ---
        left_layout = QVBoxLayout()
        
        sel_layout = QVBoxLayout()
        sel_layout.addWidget(self.epoch_label)
        sel_layout.addWidget(self.epoch_combo)
        sel_layout.addSpacing(8)
        sel_layout.addWidget(self.file_label)
        sel_layout.addWidget(self.file_combo)
        sel_layout.addSpacing(12)
        sel_layout.addWidget(self.update_button)
        sel_layout.addStretch(1)
        
        self.data_group.setLayout(sel_layout)
        left_layout.addWidget(self.data_group)
        left_layout.addStretch(1)

        # Container for left side
        left_widget = QWidget()
        left_widget.setLayout(left_layout)
        left_widget.setMaximumWidth(250)

        # --- Right Side Layout ---
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)
        
        right_widget = QWidget()
        right_widget.setLayout(right_layout)

        # --- Main Layout ---
        main_layout = QHBoxLayout()
        main_layout.addWidget(left_widget)
        main_layout.addWidget(right_widget)
        self.setLayout(main_layout)

    # ------------------------------------------------------------------ #
    #                            LOGIC                                   #
    # ------------------------------------------------------------------ #
    def update_files(self):
        year = self.epoch_combo.currentText()
        if year in ms.file_num:
            n = ms.file_num[year] + 1
            self.file_combo.clear()
            self.file_combo.addItems([str(i) for i in range(n)])

    def update_plot(self):
        year = self.epoch_combo.currentText()
        try:
            filenum = int(self.file_combo.currentText())
        except ValueError:
            return

        self.ax.clear()
        
        # Extract U/V Coordinates
        # ms.Extract returns an array of shape (6, n_wavelengths) usually, 
        # or (6,) if not spectral. Based on your snippet, we treat it as per baseline.
        try:
            UCOORD = ms.Extract(year, filenum, dataType="UCOORD")
            VCOORD = ms.Extract(year, filenum, dataType="VCOORD")
        except Exception as e:
            self.ax.text(0.5, 0.5, f"Error extracting UV data:\n{e}", ha='center')
            self.canvas.draw()
            return

        # Plot settings
        self.ax.set_title(f"U-V Coverage | {year} | File {filenum}")
        self.ax.set_xlabel("U [m]")
        self.ax.set_ylabel("V [m]")
        self.ax.grid(True, alpha=0.3, linestyle=":")
        self.ax.set_aspect('equal') # Crucial for UV maps

        # Plot loop
        for baseline in range(6):
            try:
                name = ms.baseline_name(year, baseline)
                color = self.baseline_color_list[baseline]
                
                # Get the specific U/V for this baseline
                # Note: UCOORD[baseline] might be a single float or an array depending on data structure
                u = UCOORD[baseline]
                v = VCOORD[baseline]

                # Plot (u, v)
                self.ax.scatter(
                    u, v, 
                    label=name, 
                    c=color, 
                    marker="s", 
                    s=25, 
                    edgecolors="black", 
                    linewidths=0.5
                )
                
                # Plot conjugate (-u, -v) - usually same color, no label to avoid duplicate legend
                self.ax.scatter(
                    -u, -v, 
                    c=color, 
                    marker="s", 
                    s=25,
                    edgecolors="black",
                    linewidths=0.5,
                    alpha=0.6 # Make conjugate slightly fainter if desired
                )
            except IndexError:
                continue

        # Limits and axes lines
        self.ax.set_xlim(-140, 140) # Slightly wider than 100 to fit labels if needed, or stick to 100
        self.ax.set_ylim(-140, 140)
        self.ax.axvline(x=0, c="black", lw=0.8, linestyle="--", alpha=0.5)
        self.ax.axhline(y=0, c="black", lw=0.8, linestyle="--", alpha=0.5)

        self.ax.legend(loc='upper right', fontsize='small', framealpha=0.9)
        self.canvas.draw()
        self._have_plot = True