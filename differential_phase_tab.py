# differential_phase_tab.py

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QComboBox, QPushButton, QCheckBox, QFormLayout
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

import numpy as np
import GUI_master_code as ms


class DifferentialPhaseTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # for zoom preservation
        self._have_plot = False

        # simple baseline colors (same order as baseline index 0..5)
        self.baseline_color_list = [
            "#5ac4ee", "#63e3bf", "#b95cf4",
            "#8cf55a", "#f559b7", "#f7b54c"
        ]

        self._create_widgets()
        self._create_layouts()

    # ------------------------------------------------------------------ #
    #                          WIDGET SETUP                              #
    # ------------------------------------------------------------------ #
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

        self.chk_plot_br_gamma = QCheckBox("Plot Brγ line")
        self.chk_plot_br_gamma.setChecked(False)
        self.chk_plot_br_gamma.setStyleSheet("font-weight: bold;")

        self.chk_tot = QCheckBox("Plot VISPHITOT")
        self.chk_tot.setChecked(False)
        self.chk_tot.setStyleSheet("font-weight: bold;")

        self.chk_cont = QCheckBox("Plot VISPHICONT")
        self.chk_cont.setChecked(False)
        self.chk_cont.setStyleSheet("font-weight: bold;")

        self.chk_philine = QCheckBox("PHILINE (pure line)")
        self.chk_philine.setChecked(False)
        self.chk_philine.setStyleSheet("font-weight: bold;")

        self.chk_philine_photo = QCheckBox("Apply photo absorption (PHILINE)")
        self.chk_philine_photo.setChecked(False)
        self.chk_philine_photo.setStyleSheet("font-weight: bold;")

        # α₁ / α₂ selection (with "None" semantics like NormFlux tab)
        self.alpha1_label = QLabel("α₁ spectrum:")
        self.alpha1_combo = QComboBox()
        self.alpha1_combo.addItem("None")
        for i in range(9):   # 0..8
            self.alpha1_combo.addItem(str(i))

        self.alpha2_label = QLabel("α₂ spectrum:")
        self.alpha2_combo = QComboBox()
        self.alpha2_combo.addItem("None")
        for i in range(5):   # 0..4
            self.alpha2_combo.addItem(str(i))

        # disabled unless photo-corr box is on
        self.alpha1_label.setEnabled(False)
        self.alpha1_combo.setEnabled(False)
        self.alpha2_label.setEnabled(False)
        self.alpha2_combo.setEnabled(False)

        self.chk_philine_photo.toggled.connect(self._update_photo_controls)

        # --- Matplotlib canvas + toolbar -----------------------------------
        self.figure = Figure(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Connect epoch change
        self.epoch_combo.currentTextChanged.connect(self.update_baselines)
        self.epoch_combo.currentTextChanged.connect(self.update_files)

        # Initialize for default epoch
        self.update_baselines()
        self.update_files()

    def _create_layouts(self):
        # --- Left: data selection inside group box ---
        left_selection_layout = QVBoxLayout()
        left_selection_layout.addWidget(self.epoch_label)
        left_selection_layout.addWidget(self.epoch_combo)
        left_selection_layout.addSpacing(8)

        left_selection_layout.addWidget(self.file_label)
        left_selection_layout.addWidget(self.file_combo)
        left_selection_layout.addSpacing(8)

        left_selection_layout.addWidget(self.baseline_label)
        left_selection_layout.addWidget(self.baseline_combo)
        left_selection_layout.addSpacing(8)

        left_selection_layout.addWidget(self.update_button)
        left_selection_layout.addStretch(1)

        self.data_group.setLayout(left_selection_layout)

        # --- Options layout ---
        options_layout = QVBoxLayout()
        options_layout.addWidget(self.chk_plot_br_gamma)
        options_layout.addWidget(self.chk_tot)
        options_layout.addWidget(self.chk_cont)
        options_layout.addWidget(self.chk_philine)
        options_layout.addWidget(self.chk_philine_photo)

        alphas_layout = QFormLayout()
        alphas_layout.addRow(self.alpha1_label, self.alpha1_combo)
        alphas_layout.addRow(self.alpha2_label, self.alpha2_combo)
        options_layout.addLayout(alphas_layout)

        options_layout.addStretch(1)
        self.options_group.setLayout(options_layout)

        # --- Left column: data + options ---
        left_column = QVBoxLayout()
        left_column.addWidget(self.data_group)
        left_column.addWidget(self.options_group)
        left_column.addStretch(1)

        left_widget = QWidget()
        left_widget.setLayout(left_column)

        # --- Right: toolbar + canvas ---
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)

        right_widget = QWidget()
        right_widget.setLayout(right_layout)

        # --- Main layout ---
        main_layout = QHBoxLayout()
        main_layout.addWidget(left_widget, stretch=0)
        main_layout.addWidget(right_widget, stretch=1)

        self.setLayout(main_layout)

    # ------------------------------------------------------------------ #
    #                        DATA COMBO HELPERS                          #
    # ------------------------------------------------------------------ #
    def update_baselines(self):
        year = self.epoch_combo.currentText()
        baseline_options = [ms.baseline_name(year=year, baseline=i) for i in range(6)]
        self.baseline_combo.clear()
        # first entry: "All baselines"
        self.baseline_combo.addItem("All baselines")
        self.baseline_combo.addItems(baseline_options)

    def update_files(self):
        year = self.epoch_combo.currentText()
        n_files = ms.file_num[year] + 1      # or ms.file_num[year], depending on convention
        self.file_combo.clear()
        self.file_combo.addItems([str(i) for i in range(n_files)])

    # ------------------------------------------------------------------ #
    #                        PHOTO-CORR HELPERS                          #
    # ------------------------------------------------------------------ #
    def _update_photo_controls(self, checked: bool):
        self.alpha1_label.setEnabled(checked)
        self.alpha1_combo.setEnabled(checked)
        self.alpha2_label.setEnabled(checked)
        self.alpha2_combo.setEnabled(checked)

    def _get_alpha_params_from_ui(self):
        """
        Same semantics as in NormFlux tab:

        - "None" -> alpha = 1 (neutralize), whichSpectrum = None
        - "k"    -> alpha = None, whichSpectrum = k (use Pollux spectrum k)
        """
        a1_text = self.alpha1_combo.currentText()
        if a1_text == "None":
            alpha1 = 1
            a1_ws  = None
        else:
            alpha1 = None
            a1_ws  = int(a1_text)

        a2_text = self.alpha2_combo.currentText()
        if a2_text == "None":
            alpha2 = 1
            a2_ws  = None
        else:
            alpha2 = None
            a2_ws  = int(a2_text)

        return alpha1, alpha2, a1_ws, a2_ws

    # ------------------------------------------------------------------ #
    #                             PLOTTING                               #
    # ------------------------------------------------------------------ #
    def update_plot(self):
        year    = self.epoch_combo.currentText()
        filenum = int(self.file_combo.currentText())

        lam = ms.WAVE(year, filenum)

        # --- ALWAYS preserve zoom if a plot exists -------------------------
        prev_xlim = prev_ylim = None
        if self._have_plot:
            prev_xlim = self.ax.get_xlim()
            prev_ylim = self.ax.get_ylim()

        self.ax.clear()
        self.ax.set_title("Differential Phase")
        self.ax.set_xlabel("Wavelength [m]")
        self.ax.set_ylabel("Differential phase [deg]")

        # --- check if "All baselines" or a single one ----------------------
        idx = self.baseline_combo.currentIndex()

        if idx == 0:
            # All baselines: clean multi-baseline view, only VISPHI
            for b in range(6):
                dphase = ms.VISPHI(year, filenum, b)  # already in degrees
                color = self.baseline_color_list[b % len(self.baseline_color_list)]
                label = ms.baseline_name(year=year, baseline=b)
                self.ax.plot(lam, dphase, color=color, label=label)

        else:
            # Single baseline: combo index 1..6 -> baseline index 0..5
            baseline_index = idx - 1
            baseline_name  = self.baseline_combo.currentText()

            # base differential phase (deg)
            dphase = ms.VISPHI(year, filenum, baseline_index)
            color = self.baseline_color_list[baseline_index % len(self.baseline_color_list)]
            self.ax.plot(lam, dphase, color=color, label=f"{baseline_name} VISPHI")

            # VISPHITOT (deg)
            if self.chk_tot.isChecked():
                phi_tot = ms.VISPHITOT(year, filenum, baseline_index)
                self.ax.plot(
                    lam, phi_tot,
                    color="black",
                    linestyle="--",
                    linewidth=1.0,
                    label=f"{baseline_name} VISPHITOT"
                )

            # VISPHICONT (deg)
            if self.chk_cont.isChecked():
                phi_cont = ms.VISPHICONT(year, filenum, baseline_index)
                self.ax.plot(
                    lam, phi_cont,
                    color="tab:gray",
                    linestyle="-.",
                    linewidth=1.0,
                    label=f"{baseline_name} VISPHICONT"
                )

            # --- PHILINE (pure-line phase) and photo-corrected PHILINE -----
            mask = ms.mask_brg(year, filenum)
            wave_line = ms.WAVE(year, filenum)[mask]

            philine_err_deg = None  # cache for reuse

            # PHILINE without photospheric correction (radians -> degrees)
            if self.chk_philine.isChecked():
                philine_rad = ms.PHILINE(
                    year=year,
                    filenum=filenum,
                    baseline=baseline_index,
                    PhotosphericCorr=False,
                    alpha1=None,
                    alpha2=None,
                    beta1=0,
                    beta2=0,
                    returnError=False
                )[mask]

                philine_err_rad = ms.PHILINE(
                    year=year,
                    filenum=filenum,
                    baseline=baseline_index,
                    PhotosphericCorr=False,
                    alpha1=None,
                    alpha2=None,
                    beta1=0,
                    beta2=0,
                    returnError=True
                )[mask]

                philine_deg = np.rad2deg(philine_rad)
                philine_err_deg = np.rad2deg(philine_err_rad)

                self.ax.errorbar(
                    wave_line,
                    philine_deg,
                    yerr=philine_err_deg,
                    fmt="s",
                    capsize=4,
                    markeredgecolor="black",
                    markerfacecolor=color,
                    ecolor=color,
                    label=f"PHILINE ({baseline_name})"
                )

            # PHILINE with photospheric correction (also radians -> degrees)
            if self.chk_philine_photo.isChecked():
                alpha1_val, alpha2_val, a1_ws, a2_ws = self._get_alpha_params_from_ui()

                philine_corr_rad = ms.PHILINE(
                    year=year,
                    filenum=filenum,
                    baseline=baseline_index,
                    PhotosphericCorr=True,
                    alpha1_whichSpectrum=a1_ws,
                    alpha2_whichSpectrum=a2_ws,
                    alpha1=alpha1_val,
                    alpha2=alpha2_val,
                    beta1=0,
                    beta2=0,
                    returnError=False
                )[mask]

                philine_corr_deg = np.rad2deg(philine_corr_rad)

                # errors (reuse if already computed above)
                if philine_err_deg is None:
                    philine_err_rad = ms.PHILINE(
                        year=year,
                        filenum=filenum,
                        baseline=baseline_index,
                        PhotosphericCorr=True,
                        alpha1_whichSpectrum=a1_ws,
                        alpha2_whichSpectrum=a2_ws,
                        alpha1=alpha1_val,
                        alpha2=alpha2_val,
                        beta1=0,
                        beta2=0,
                        returnError=True
                    )[mask]
                    philine_err_deg = np.rad2deg(philine_err_rad)

                label_str = "PHILINE + photo abs ("
                if a1_ws is None:
                    label_str += "α₁=None"
                else:
                    label_str += f"α₁={a1_ws}"
                label_str += ", "
                if a2_ws is None:
                    label_str += "α₂=None"
                else:
                    label_str += f"α₂={a2_ws}"
                label_str += ")"

                self.ax.errorbar(
                    wave_line,
                    philine_corr_deg,
                    yerr=philine_err_deg,
                    fmt="s",
                    capsize=4,
                    markeredgecolor="black",
                    markerfacecolor="none",
                    ecolor=color,
                    label=label_str
                )

        # --- optional Brγ line --------------------------------------------
        if self.chk_plot_br_gamma.isChecked():
            self.ax.axvline(
                x=ms.bg,
                color="tab:red",
                linestyle="--",
                linewidth=1.2,
                label="Brγ (2.16612 μm)"
            )

        # --- restore zoom --------------------------------------------------
        if prev_xlim is not None and prev_ylim is not None:
            self.ax.set_xlim(prev_xlim)
            self.ax.set_ylim(prev_ylim)

        self.ax.legend()
        self.canvas.draw()

        self._have_plot = True
