import numpy as np
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QComboBox, QPushButton, QCheckBox, QFormLayout, QPlainTextEdit
)
from PyQt5.QtGui import QFont

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

import GUI_master_code as ms


class NormalizedFluxTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # for zoom preservation
        self._last_view_params = None   # (year, filenum, baseline_index)
        self._have_plot = False

        self._create_widgets()
        self._create_layouts()

    def _create_widgets(self):
        # --- Left side: Data selection group -------------------------------
        self.data_group = QGroupBox("Data selection")

        # Epoch
        self.epoch_label = QLabel("Epoch:")
        self.epoch_label.setStyleSheet("font-weight: bold;")
        self.epoch_combo = QComboBox()
        self.epoch_combo.addItems([
            "2018_red", "2019_red", "2020_red",
            "2021_red", "2022a_red", "2022b_red", "2023_red"
        ])

        # File number
        self.file_label = QLabel("File number:")
        self.file_label.setStyleSheet("font-weight: bold;")
        self.file_combo = QComboBox()

        # Baseline (not used in NormFlux physics, but present for consistency)
        self.baseline_label = QLabel("Baseline:")
        self.baseline_label.setStyleSheet("font-weight: bold;")
        self.baseline_combo = QComboBox()

        # Update button
        self.update_button = QPushButton("Update plot")
        self.update_button.clicked.connect(self.update_plot)

        # --- Options group -------------------------------------------------
        self.options_group = QGroupBox("Options")

        # 0) Brγ vertical line
        self.chk_brg = QCheckBox("Plot Brγ line")
        self.chk_brg.setChecked(False)
        self.chk_brg.setStyleSheet("font-weight: bold;")

        # 1) Tellurics
        self.chk_tellurics = QCheckBox("Plot tellurics")
        self.chk_tellurics.setChecked(False)
        self.chk_tellurics.setStyleSheet("font-weight: bold;")

        # 2) FLCcorr (here: NormFlux(returnFLC=True) + error band)
        self.chk_flc_corr = QCheckBox("Plot FLCcorr (+ error)")
        self.chk_flc_corr.setChecked(False)
        self.chk_flc_corr.setStyleSheet("font-weight: bold;")

        # 3) Photo-corrected flux (using alpha1/alpha2 choices)
        self.chk_photo_corr = QCheckBox("Plot photo-corrected flux")
        self.chk_photo_corr.setChecked(False)
        self.chk_photo_corr.setStyleSheet("font-weight: bold;")

        # α₁ / α₂ selection with "None" option
        self.alpha1_label = QLabel("α₁ spectrum:")
        self.alpha1_combo = QComboBox()
        self.alpha1_combo.addItem("None")
        for i in range(9):                 # 0..8
            self.alpha1_combo.addItem(str(i))

        self.alpha2_label = QLabel("α₂ spectrum:")
        self.alpha2_combo = QComboBox()
        self.alpha2_combo.addItem("None")
        for i in range(5):                 # 0..4
            self.alpha2_combo.addItem(str(i))

        # disable α-controls unless photo-corr is checked
        self.alpha1_label.setEnabled(False)
        self.alpha1_combo.setEnabled(False)
        self.alpha2_label.setEnabled(False)
        self.alpha2_combo.setEnabled(False)

        self.chk_photo_corr.toggled.connect(self._update_photo_controls)

        # 4) Absorption 1 and 2
        self.chk_abs1 = QCheckBox("Plot absorption 1")
        self.chk_abs1.setChecked(False)
        self.chk_abs1.setStyleSheet("font-weight: bold;")

        self.chk_abs2 = QCheckBox("Plot absorption 2")
        self.chk_abs2.setChecked(False)
        self.chk_abs2.setStyleSheet("font-weight: bold;")

        # --- Gaussian fit controls -----------------------------------------
        self.gauss_fit_button = QPushButton("Fit double Gaussian (Brγ)")
        self.gauss_fit_button.setStyleSheet("font-weight: bold;")
        self.gauss_fit_button.clicked.connect(self.run_gaussian_fit)

        # --- Terminal-like output for fit results --------------------------
        self.fit_output = QPlainTextEdit()
        self.fit_output.setReadOnly(True)
        self.fit_output.setObjectName("FitOutput")
        self.fit_output.setMinimumHeight(150)
        self.fit_output.setLineWrapMode(QPlainTextEdit.NoWrap)

        font = QFont("Courier New")
        font.setPointSize(9)
        self.fit_output.setFont(font)
        self.fit_output.setPlainText(
            "Fit report will appear here.\n"
            "Click 'Fit double Gaussian (Brγ)' to run a fit."
        )

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
        # --- Data selection layout (inside data_group) ---
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

        # --- Options layout (inside options_group) ---
        options_layout = QVBoxLayout()
        options_layout.addWidget(self.chk_brg)
        options_layout.addWidget(self.chk_tellurics)
        options_layout.addWidget(self.chk_flc_corr)
        options_layout.addWidget(self.chk_photo_corr)

        alphas_layout = QFormLayout()
        alphas_layout.addRow(self.alpha1_label, self.alpha1_combo)
        alphas_layout.addRow(self.alpha2_label, self.alpha2_combo)
        options_layout.addLayout(alphas_layout)

        options_layout.addWidget(self.chk_abs1)
        options_layout.addWidget(self.chk_abs2)

        # Gaussian fit button
        options_layout.addWidget(self.gauss_fit_button)

        options_layout.addStretch(1)
        self.options_group.setLayout(options_layout)

        # --- Left column: data group + options group ---
        left_column = QVBoxLayout()
        left_column.addWidget(self.data_group)
        left_column.addWidget(self.options_group)
        left_column.addStretch(1)

        left_widget = QWidget()
        left_widget.setLayout(left_column)

        # --- Right side: toolbar + canvas + fit output ---
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas, stretch=1)
        right_layout.addWidget(self.fit_output, stretch=0)

        right_widget = QWidget()
        right_widget.setLayout(right_layout)

        # --- Main layout: left and right ---
        main_layout = QHBoxLayout()
        main_layout.addWidget(left_widget, stretch=0)
        main_layout.addWidget(right_widget, stretch=1)

        self.setLayout(main_layout)

    # --- helpers -----------------------------------------------------------
    def _update_photo_controls(self, checked: bool):
        self.alpha1_label.setEnabled(checked)
        self.alpha1_combo.setEnabled(checked)
        self.alpha2_label.setEnabled(checked)
        self.alpha2_combo.setEnabled(checked)

    def _get_alpha_params_from_ui(self):
        """
        Interpret the alpha1/alpha2 combo boxes in your semantics:

        - "None" -> alpha = 1 (neutralize this component), whichSpectrum = None
        - "k"    -> alpha = None, whichSpectrum = k (use Pollux spectrum k)
        """
        # α1
        a1_text = self.alpha1_combo.currentText()
        if a1_text == "None":
            alpha1 = 1
            a1_ws = None
        else:
            alpha1 = None
            a1_ws = int(a1_text)

        # α2
        a2_text = self.alpha2_combo.currentText()
        if a2_text == "None":
            alpha2 = 1
            a2_ws = None
        else:
            alpha2 = None
            a2_ws = int(a2_text)

        return alpha1, alpha2, a1_ws, a2_ws

    # --- shared logic: similar to visibility tab --------------------------
    def update_baselines(self):
        year = self.epoch_combo.currentText()
        baseline_options = [ms.baseline_name(year=year, baseline=i) for i in range(6)]
        self.baseline_combo.clear()
        self.baseline_combo.addItems(baseline_options)

    def update_files(self):
        year = self.epoch_combo.currentText()
        n_files = ms.file_num[year] + 1
        self.file_combo.clear()
        self.file_combo.addItems([str(i) for i in range(n_files)])

    def run_gaussian_fit(self):
        """Run double-Gaussian fit and print a scrollable report, then overplot in green."""
        year = self.epoch_combo.currentText()
        try:
            filenum = int(self.file_combo.currentText())
        except ValueError:
            self.fit_output.setPlainText("No valid file selected.")
            return

        # execute fit without plotting or Rich console spam
        try:
            res = ms.FitGaussian(
                year=year,
                filenum=filenum,
                showHelp=False,
                printFitParameters=False,
                plotLine=None,
                do_plot=False,
                return_dict=True,
            )
        except Exception as e:
            self.fit_output.setPlainText(f"Fit failed:\n{repr(e)}")
            return

        # build a compact, terminal-like report
        lines = []
        lines.append(f"Double-Gaussian fit   epoch={year}   file={filenum}")
        lines.append("-" * 60)
        lines.append("Fit parameters (Å):")
        lines.append(f"  A1   = {res['A1']:.4g} ± {res['A1_err']:.2g}")
        lines.append(f"  μ1   = {res['mu1_A']:.4f} ± {res['mu1_A_err']:.4f} Å")
        lines.append(f"  σ1   = {res['s1_A']:.4f} ± {res['s1_A_err']:.4f} Å")
        lines.append(f"  A2   = {res['A2']:.4g} ± {res['A2_err']:.2g}")
        lines.append(f"  μ2   = {res['mu2_A']:.4f} ± {res['mu2_A_err']:.4f} Å")
        lines.append(f"  σ2   = {res['s2_A']:.4f} ± {res['s2_A_err']:.4f} Å")
        lines.append(f"  C    = {res['C']:.5f} ± {res['C_err']:.2g}")
        lines.append("")
        lines.append("Derived (Å, m, km/s):")
        lines.append(f"  FWHM1  = {res['FWHM1_A']:.4f} ± {res['FWHM1_A_err']:.4f} Å")
        lines.append(f"  FWHM2  = {res['FWHM2_A']:.4f} ± {res['FWHM2_A_err']:.4f} Å")
        lines.append(f"  Δμ     = {res['dmu_A']:.4f} ± {res['dmu_A_err']:.4f} Å")
        lines.append(f"  μ̄      = {res['mid_A']:.4f} ± {res['mid_A_err']:.4f} Å")
        lines.append("")
        lines.append(f"  Δv     = {res['delta_v']:.2f} ± {res['delta_v_err']:.2f} km/s")
        lines.append(f"  v_mid  = {res['v_mid']:.2f} ± {res['v_mid_err']:.2f} km/s")
        lines.append("")
        lines.append("Equivalent widths:")
        lines.append(f"  EW₁    = {res['EW1_A']:.3f} ± {res['EW1_A_err']:.2f} Å")
        lines.append(f"  EW₂    = {res['EW2_A']:.3f} ± {res['EW2_A_err']:.2f} Å")
        lines.append(f"  EW_tot = {res['EW_total_A']:.3f} ± {res['EW_total_A_err']:.2f} Å")
        lines.append("")
        lines.append(f"V/R     = {res['V_over_R']:.3f} ± {res['V_over_R_err']:.2f}")
        lines.append(f"χ²_red  = {res['chi2_red']:.3f}")

        text = "\n".join(lines)
        self.fit_output.setPlainText(text)
        self.fit_output.verticalScrollBar().setValue(0)

        # Overplot the fitted double Gaussian on the panel plot (in green)
        try:
            WL = ms.WAVE(year, filenum)
            mask = (WL >= 2.161e-6) & (WL <= 2.172e-6)

            xx_m = np.linspace(WL[mask].min(), WL[mask].max(), 1500)
            xx_A = xx_m * 1e10  # convert to Å

            def double_gaussian_with_offset(x, A1, mu1, s1, A2, mu2, s2, C):
                g1 = A1 * np.exp(-(x - mu1) ** 2 / (2 * s1 ** 2))
                g2 = A2 * np.exp(-(x - mu2) ** 2 / (2 * s2 ** 2))
                return g1 + g2 + C

            yy_fit = double_gaussian_with_offset(
                xx_A,
                res["A1"], res["mu1_A"], res["s1_A"],
                res["A2"], res["mu2_A"], res["s2_A"],
                res["C"]
            )

            for line in list(self.ax.lines):
                if line.get_label() == "gauss_fit":
                    line.remove()

            self.ax.plot(
                xx_m, yy_fit,
                color="green",
                lw=1.2,
                label="gauss_fit"
            )

            handles, labels = self.ax.get_legend_handles_labels()
            self.ax.legend(handles, labels, frameon=False)

            self.canvas.draw()

        except Exception as e:
            self.fit_output.appendPlainText(f"\nPlot overlay failed: {repr(e)}")

    def update_plot(self):
        year = self.epoch_combo.currentText()
        filenum = int(self.file_combo.currentText())
        baseline_index = self.baseline_combo.currentIndex()  # not used in physics

        wl = ms.WAVE(year, filenum)

        # --- Pollux-corrected wavelength for absorption spectra -------------
        # Uses the low-resolution wavelength grid shifted to align the Pollux
        # absorption with the Brγ feature, as defined in PolluxSpectra.
        try:
            wave_pollux = ms.PolluxSpectra(
                year,
                filenum,
                showHelp=False,
                whichStar="Primary",
                whichSpectrum=0,
                fullResSpectra=True,
                plot=False,
                listSpectra=False,
                returnSpectra=True,
                returnLowResWave=True,
            )
        except Exception:
            # Fallback if Pollux files are missing or fail to load
            wave_pollux = wl

        # --- zoom preservation logic --------------------------------------
        current_params = (year, filenum, baseline_index)
        preserve_view = (self._have_plot and
                         self._last_view_params == current_params)

        prev_xlim = prev_ylim = None
        if preserve_view:
            prev_xlim = self.ax.get_xlim()
            prev_ylim = self.ax.get_ylim()

        # --- main raw normalized flux (RAW_NOWM_FLUX) ----------------------
        flux_raw = ms.NormFlux(
            year=year,
            filenum=filenum,
            PhotosphericCorr=False,
            alpha1=None,
            alpha2=None,
            alpha1_whichSpectrum=None,
            alpha2_whichSpectrum=None,
            returnFLC=False,
            returnFLCcorr=False,
            returnPhotoCorrectionFlux=False,
            returnRAWFLUX=True,
            returnTellurics=False,
            returnAbsorbtion1=False,
            returnAbsorbtion2=False
        )

        self.ax.clear()
        self.ax.set_title("Normalized Flux")
        self.ax.set_xlabel("Wavelength [m]")
        self.ax.set_ylabel("Normalized flux")

        self.ax.plot(
            wl,
            flux_raw,
            color="orange",
            linewidth=1.5,
            label="Raw normalized flux"
        )

        # --- optional Brγ vertical line -----------------------------------
        if self.chk_brg.isChecked():
            self.ax.axvline(
                x=ms.bg,
                color="tab:red",
                linestyle="--",
                linewidth=1.2,
                label="Brγ (2.16612 μm)"
            )

        # --- optional FLCcorr + error band --------------------------------
        if self.chk_flc_corr.isChecked():
            normFlux = ms.NormFlux(year, filenum, returnFLC=True)
            normFluxErr = ms.NormFlux(year, filenum, returnFLCerr=True)
            wave = wl

            self.ax.fill_between(
                wave,
                normFlux - normFluxErr,
                normFlux + normFluxErr,
                alpha=0.2,
                color="black",
                label="FLC ± err"
            )

            self.ax.plot(
                wave,
                normFlux,
                color="black",
                linewidth=1.2,
                label="FLC"
            )

        # --- alpha parameters from UI (used below) ------------------------
        alpha1_val, alpha2_val, a1_ws, a2_ws = self._get_alpha_params_from_ui()

        # --- optional photo-corrected flux --------------------------------
        if self.chk_photo_corr.isChecked():
            flux_photo = ms.NormFlux(
                year=year,
                filenum=filenum,
                PhotosphericCorr=True,
                alpha1_whichSpectrum=a1_ws,
                alpha2_whichSpectrum=a2_ws,
                alpha1=alpha1_val,
                alpha2=alpha2_val,
                returnFLC=False,
                returnFLCcorr=False,
                returnPhotoCorrectionFlux=True,
                returnRAWFLUX=False,
                returnTellurics=False,
                returnAbsorbtion1=False,
                returnAbsorbtion2=False
            )

            label_str = "Photo-corr ("
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

            # Keep photo-corrected flux on the original wavelength grid
            self.ax.plot(
                wl, flux_photo,
                color="tab:red",
                linestyle="--",
                linewidth=1,
                label=label_str
            )

        # --- optional tellurics overlay -----------------------------------
        if self.chk_tellurics.isChecked():
            tell = ms.NormFlux(
                year=year,
                filenum=filenum,
                PhotosphericCorr=False,
                alpha1=None,
                alpha2=None,
                alpha1_whichSpectrum=None,
                alpha2_whichSpectrum=None,
                returnTellurics=True,
                returnRAWFLUX=False,
                returnFLC=False,
                returnFLCcorr=False,
                returnPhotoCorrectionFlux=False,
                returnAbsorbtion1=False,
                returnAbsorbtion2=False
            )
            self.ax.plot(
                wl, tell,
                color="blue",
                alpha=0.7,
                label="Tellurics"
            )

        # --- optional absorption features (Pollux-based) -------------------
        # These are the Pollux absorption spectra; plot them on the
        # Pollux-corrected wavelength axis `wave_pollux`.
        if self.chk_abs1.isChecked():
            abs1 = ms.NormFlux(
                year=year,
                filenum=filenum,
                PhotosphericCorr=True,
                alpha1_whichSpectrum=a1_ws,
                alpha2_whichSpectrum=a2_ws,
                alpha1=alpha1_val,
                alpha2=alpha2_val,
                returnFLC=False,
                returnFLCcorr=False,
                returnPhotoCorrectionFlux=False,
                returnRAWFLUX=False,
                returnTellurics=False,
                returnAbsorbtion1=True,
                returnAbsorbtion2=False
            )
            self.ax.plot(
                wave_pollux, abs1,
                color="tab:purple",
                linestyle="-.",
                linewidth=1,
                label="Absorption 1"
            )

        if self.chk_abs2.isChecked():
            abs2 = ms.NormFlux(
                year=year,
                filenum=filenum,
                PhotosphericCorr=True,
                alpha1_whichSpectrum=a1_ws,
                alpha2_whichSpectrum=a2_ws,
                alpha1=alpha1_val,
                alpha2=alpha2_val,
                returnFLC=False,
                returnFLCcorr=False,
                returnPhotoCorrectionFlux=False,
                returnRAWFLUX=False,
                returnTellurics=False,
                returnAbsorbtion1=False,
                returnAbsorbtion2=True
            )
            self.ax.plot(
                wave_pollux, abs2,
                color="tab:orange",
                linestyle="-.",
                linewidth=1,
                label="Absorption 2"
            )

        # --- default x/y limits if first time ------------------------------
        if not preserve_view:
            self.ax.set_xlim(2.16e-6, 2.172e-6)
            self.ax.set_ylim(0.7, 1.4)

        # --- restore zoom if same dataset ---------------------------------
        if preserve_view and prev_xlim is not None and prev_ylim is not None:
            self.ax.set_xlim(prev_xlim)
            self.ax.set_ylim(prev_ylim)

        self.ax.legend(frameon=False)
        self.canvas.draw()

        self._last_view_params = current_params
        self._have_plot = True
