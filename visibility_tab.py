import sys
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QComboBox, QPushButton, QCheckBox, QFormLayout, QToolButton
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QColorDialog
from PyQt5.QtCore import Qt

import matplotlib.pyplot as plt
import numpy as np
import GUI_master_code as ms


class VisibilityTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # for zoom preservation
        self._have_plot = False

        # list of baseline colors (order must match your baseline indices 0..5)
        self.baseline_color_list = [
            "#5ac4ee", "#63e3bf", "#b95cf4",
            "#8cf55a", "#f559b7", "#f7b54c"
        ]

        # color lock flags
        self.visline_color_locked = False
        self.photoabs_color_locked = False

        self._create_widgets()
        self._create_layouts()

    def _create_widgets(self):
        # --- Left side: Data selection group -------------------------------
        self.data_group = QGroupBox("Data selection")

        # Epoch selector
        self.epoch_label = QLabel("Epoch:")
        self.epoch_label.setStyleSheet("font-weight: bold;")
        self.epoch_combo = QComboBox()
        self.epoch_combo.addItems([
            "2018_red", "2019_red", "2020_red",
            "2021_red", "2022a_red", "2022b_red", "2023_red"
        ])

        # File selector (combo)
        self.file_label = QLabel("File number:")
        self.file_label.setStyleSheet("font-weight: bold;")
        self.file_combo = QComboBox()

        # Baseline selector
        self.baseline_label = QLabel("Baseline:")
        self.baseline_label.setStyleSheet("font-weight: bold;")
        self.baseline_combo = QComboBox()

        # Button
        self.update_button = QPushButton("Update plot")
        self.update_button.clicked.connect(self.update_plot)

        # ---- Options group ----
        self.options_group = QGroupBox("Options")

        self.chk_plot_br_gamma = QCheckBox("Plot Bracket γ")
        self.chk_plot_br_gamma.setChecked(False)
        self.chk_plot_br_gamma.setStyleSheet("font-weight: bold;")

        self.chk_continuum = QCheckBox("Plot continuum")
        self.chk_continuum.setChecked(False)
        self.chk_continuum.setStyleSheet("font-weight: bold;")

        self.chk_visline = QCheckBox("VISLINE")
        self.chk_visline.setChecked(False)
        self.chk_visline.setStyleSheet("font-weight: bold;")

        self.chk_photoabs = QCheckBox("Apply photo absorption")
        self.chk_photoabs.setChecked(False)
        self.chk_photoabs.setStyleSheet("font-weight: bold;")

        # --- Color state (defaults based on baseline 0) ---
        default_hex = self.baseline_color_list[0]
        self.visline_color = QColor(default_hex)
        self.photoabs_color = QColor(default_hex)

        # small color buttons
        self.btn_visline_color = QToolButton()
        self.btn_visline_color.setToolTip("Choose VISLINE color")
        self._update_button_color(self.btn_visline_color, self.visline_color)
        self.btn_visline_color.clicked.connect(self.choose_visline_color)

        self.btn_photoabs_color = QToolButton()
        self.btn_photoabs_color.setToolTip("Choose photo absorption color")
        self._update_button_color(self.btn_photoabs_color, self.photoabs_color)
        self.btn_photoabs_color.clicked.connect(self.choose_photoabs_color)

        # sub-options for photo absorption
        self.alpha1_label = QLabel("α₁ spectrum:")
        self.alpha1_combo = QComboBox()
        self.alpha1_combo.addItems([str(i) for i in range(9)])  # 0..8

        self.alpha2_label = QLabel("α₂ spectrum:")
        self.alpha2_combo = QComboBox()
        self.alpha2_combo.addItems([str(i) for i in range(5)])  # 0..4

        # start disabled; enabled only when checkbox is on
        self.alpha1_label.setEnabled(False)
        self.alpha1_combo.setEnabled(False)
        self.alpha2_label.setEnabled(False)
        self.alpha2_combo.setEnabled(False)

        # connect toggle to enabling/disabling suboptions
        self.chk_photoabs.toggled.connect(self._update_photoabs_controls)

        # --- Matplotlib canvas + toolbar -----------------------------------
        self.figure = Figure(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Connect changes AFTER combos exist
        self.epoch_combo.currentTextChanged.connect(self.update_baselines)
        self.epoch_combo.currentTextChanged.connect(self.update_files)
        self.baseline_combo.currentIndexChanged.connect(self._update_default_colors)

        # Initialize for default epoch (first item)
        self.update_baselines()
        self.update_files()
        self._update_default_colors()

    def _create_layouts(self):
        # --- Left side: controls inside the "Data selection" group box ---
        left_layout = QVBoxLayout()
        left_layout.addWidget(self.epoch_label)
        left_layout.addWidget(self.epoch_combo)
        left_layout.addSpacing(8)

        left_layout.addWidget(self.file_label)
        left_layout.addWidget(self.file_combo)
        left_layout.addSpacing(8)

        left_layout.addWidget(self.baseline_label)
        left_layout.addWidget(self.baseline_combo)
        left_layout.addSpacing(8)

        left_layout.addWidget(self.update_button)
        left_layout.addStretch(1)

        self.data_group.setLayout(left_layout)

        # --- Options layout ---
        options_layout = QVBoxLayout()
        options_layout.addWidget(self.chk_plot_br_gamma)
        options_layout.addWidget(self.chk_continuum)

        # VISLINE row: checkbox + color button
        visline_row = QHBoxLayout()
        visline_row.addWidget(self.chk_visline)
        visline_row.addWidget(self.btn_visline_color)
        visline_row.addStretch(1)
        options_layout.addLayout(visline_row)

        # Photo absorption row: checkbox + color button
        photoabs_row = QHBoxLayout()
        photoabs_row.addWidget(self.chk_photoabs)
        photoabs_row.addWidget(self.btn_photoabs_color)
        photoabs_row.addStretch(1)
        options_layout.addLayout(photoabs_row)

        # alpha controls below
        alphas_layout = QFormLayout()
        alphas_layout.addRow(self.alpha1_label, self.alpha1_combo)
        alphas_layout.addRow(self.alpha2_label, self.alpha2_combo)

        options_layout.addLayout(alphas_layout)
        options_layout.addStretch(1)
        self.options_group.setLayout(options_layout)

        # --- Right side: toolbar + canvas ---
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)

        right_widget = QWidget()
        right_widget.setLayout(right_layout)

        # --- Main layout: left group + right widget ---
        main_layout = QHBoxLayout()

        left_column = QVBoxLayout()
        left_column.addWidget(self.data_group)
        left_column.addWidget(self.options_group)
        left_column.addStretch(1)

        left_widget = QWidget()
        left_widget.setLayout(left_column)
        main_layout.addWidget(left_widget, stretch=0)

        main_layout.addWidget(right_widget, stretch=1)

        self.setLayout(main_layout)

    # ------------------------------------------------------------------ #
    #                    COLOR / OPTION HELPERS                          #
    # ------------------------------------------------------------------ #
    def _update_default_colors(self):
        """
        When baseline selection changes, update default colors for VISLINE
        and photoabs *only* if the user did not pick custom colors.
        """
        idx = self.baseline_combo.currentIndex()
        # index 0 = "All baselines" (we default to baseline 0 color)
        baseline_index = 0 if idx <= 0 else idx - 1
        hexcolor = self.baseline_color_list[baseline_index]
        qc = QColor(hexcolor)

        if not self.visline_color_locked:
            self.visline_color = qc
            self._update_button_color(self.btn_visline_color, qc)

        if not self.photoabs_color_locked:
            self.photoabs_color = qc
            self._update_button_color(self.btn_photoabs_color, qc)

    def _update_button_color(self, button: QToolButton, color: QColor):
        r, g, b, _ = color.getRgb()
        button.setStyleSheet(
            f"background-color: rgb({r},{g},{b}); "
            "border: 1px solid #444; width: 20px; height: 20px;"
        )

    def choose_visline_color(self):
        color = QColorDialog.getColor(self.visline_color, self, "Choose VISLINE color")
        if color.isValid():
            self.visline_color = color
            self.visline_color_locked = True
            self._update_button_color(self.btn_visline_color, color)
            self.update_plot()

    def choose_photoabs_color(self):
        color = QColorDialog.getColor(self.photoabs_color, self, "Choose photo absorption color")
        if color.isValid():
            self.photoabs_color = color
            self.photoabs_color_locked = True
            self._update_button_color(self.btn_photoabs_color, color)
            self.update_plot()

    def _update_photoabs_controls(self, checked: bool):
        self.alpha1_label.setEnabled(checked)
        self.alpha1_combo.setEnabled(checked)
        self.alpha2_label.setEnabled(checked)
        self.alpha2_combo.setEnabled(checked)

    # ------------------------------------------------------------------ #
    #                      DATA COMBO UPDATES                            #
    # ------------------------------------------------------------------ #
    def update_baselines(self):
        year = self.epoch_combo.currentText()
        baseline_options = [ms.baseline_name(year=year, baseline=i) for i in range(6)]
        self.baseline_combo.clear()
        self.baseline_combo.addItem("All baselines")
        self.baseline_combo.addItems(baseline_options)

    def update_files(self):
        year = self.epoch_combo.currentText()
        n_files = ms.file_num[year] + 1    # keep your convention
        self.file_combo.clear()
        self.file_combo.addItem("All files")
        self.file_combo.addItems([str(i) for i in range(n_files)])

    # ------------------------------------------------------------------ #
    #                           PLOTTING                                 #
    # ------------------------------------------------------------------ #
    def update_plot(self):
        year = self.epoch_combo.currentText()

        # preserve zoom ALWAYS if we already plotted something
        prev_xlim = prev_ylim = None
        if self._have_plot:
            prev_xlim = self.ax.get_xlim()
            prev_ylim = self.ax.get_ylim()

        # determine which files to plot
        file_idx = self.file_combo.currentIndex()
        if file_idx == 0:
            # All files
            n_files = ms.file_num[year] + 1
            file_indices = list(range(n_files))
        else:
            file_indices = [int(self.file_combo.currentText())]

        # determine which baselines to plot
        base_idx = self.baseline_combo.currentIndex()
        if base_idx == 0:
            # All baselines
            baseline_indices = list(range(6))
        else:
            baseline_indices = [base_idx - 1]

        # single/single mode: full feature set; otherwise: overview only
        single_file = (len(file_indices) == 1)
        single_baseline = (len(baseline_indices) == 1)
        detailed_mode = single_file and single_baseline

        self.ax.clear()
        self.ax.set_title("Visibility amplitude")
        self.ax.set_xlabel("Wavelength [m]")
        self.ax.set_ylabel("Visibility amplitude")

        # main plotting loops
        for filenum in file_indices:
            lam = ms.WAVE(year, filenum)

            for baseline in baseline_indices:
                visamp = ms.VISAMP(year, filenum, baseline=baseline)

                color = self.baseline_color_list[baseline % len(self.baseline_color_list)]
                baseline_name = ms.baseline_name(year=year, baseline=baseline)

                if len(file_indices) > 1:
                    label = f"{baseline_name}, file {filenum}"
                else:
                    label = baseline_name

                self.ax.plot(lam, visamp, color=color, label=label)

        # Only in detailed mode (one file, one baseline) apply the extra layers
        if detailed_mode:
            filenum = file_indices[0]
            baseline = baseline_indices[0]
            baseline_name = ms.baseline_name(year=year, baseline=baseline)

            # colors for points
            visline_color_rgb = tuple(self.visline_color.getRgbF()[:3])
            photoabs_color_rgb = tuple(self.photoabs_color.getRgbF()[:3])

            # continuum
            if self.chk_continuum.isChecked():
                viscont = ms.VISAMPCONT(year, filenum, baseline)
                wl_cont = ms.WAVE(year, filenum)
                self.ax.plot(
                    wl_cont, viscont,
                    linestyle="--",
                    color="black",
                    label=f"{baseline_name} continuum"
                )

            # Brγ line
            if self.chk_plot_br_gamma.isChecked():
                self.ax.axvline(
                    x=ms.bg,
                    label="Brγ (2.16612 μm)",
                    color="red"
                )

            # mask region and wave for line-related things
            mask = ms.mask_brg(year, filenum)
            wave = ms.WAVE(year, filenum)[mask]

            # --- VISLINE (no photo absorption) ---
            visline_err = None
            if self.chk_visline.isChecked():
                visline = ms.VISLINE(year, filenum, baseline)[mask]
                visline_err = ms.VISLINE(
                    year,
                    filenum,
                    baseline,
                    returnError=True
                )[mask]

                self.ax.errorbar(
                    wave,
                    visline,
                    yerr=visline_err,
                    fmt="s",
                    capsize=4,
                    markeredgecolor="black",
                    markerfacecolor=visline_color_rgb,
                    ecolor=visline_color_rgb,
                    label=f"VISLINE ({baseline_name})"
                )

            # --- VISLINE with photo absorption ---
            if self.chk_photoabs.isChecked():
                alpha1 = int(self.alpha1_combo.currentText())
                alpha2 = int(self.alpha2_combo.currentText())

                vis_line_pa = ms.VISLINE(
                    year,
                    filenum,
                    baseline,
                    alpha1_whichSpectrum=alpha1,
                    alpha2_whichSpectrum=alpha2,
                    PhotosphericCorr=True,
                )[mask]

                if visline_err is None:
                    visline_err = ms.VISLINE(
                        year,
                        filenum,
                        baseline,
                        returnError=True
                    )[mask]

                self.ax.errorbar(
                    wave,
                    vis_line_pa,
                    yerr=visline_err,
                    fmt="s",
                    capsize=4,
                    markeredgecolor="black",
                    markerfacecolor=photoabs_color_rgb,
                    ecolor=photoabs_color_rgb,
                    label=f"VISLINE + photo abs (α₁={alpha1}, α₂={alpha2})"
                )

        # restore zoom if we had a previous view
        if prev_xlim is not None and prev_ylim is not None:
            self.ax.set_xlim(prev_xlim)
            self.ax.set_ylim(prev_ylim)
        else:
            self.ax.set_ylim(0, 1.4)

        self.ax.legend()
        self.canvas.draw()

        self._have_plot = True
