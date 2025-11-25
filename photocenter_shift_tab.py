# photocenter_shift_tab.py

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QComboBox, QPushButton, QCheckBox, QFormLayout, QSpinBox, QSlider, QFrame
)
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit  # Required for the new plot fitting
import rebound  # Required for local orbit calculation

import GUI_master_code as ms


class DoubleClickResetSlider(QSlider):
    """
    A specific slider that resets to 0 when double-clicked.
    """
    def __init__(self, orientation, parent=None):
        super().__init__(orientation, parent)

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.setValue(0)
            event.accept()
        else:
            super().mouseDoubleClickEvent(event)


class PhotocenterShiftTab(QWidget):
    """
    Tab for pure-line photocenter shift analysis.

    Modes:
    1. 2D Map (Δx, Δy): Photocenter shift in RA/Dec plane colored by velocity.
    2. Projected Shift vs PA: Pure-line photocenter shift projected onto baseline vs Position Angle.
    3. Side-by-Side: Both plots displayed simultaneously.
    """
    def __init__(self, parent=None):
        super().__init__(parent)

        self._have_plot = False
        self._len_v = None
        self._cbar = None

        # Scale for FOV sliders: value 100 -> 1.00 mas
        self._slider_scale = 100.0

        self._create_widgets()
        self._create_layouts()

    # ------------------------------------------------------------------ #
    #                          WIDGET SETUP                              #
    # ------------------------------------------------------------------ #
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

        self.baseline_label = QLabel("Baseline (Fit Set):")
        self.baseline_combo = QComboBox()

        self.update_button = QPushButton("Update Plot")
        self.update_button.setStyleSheet("font-weight: bold; padding: 5px;")
        self.update_button.clicked.connect(self.update_plot)

        # --- View Mode Group ----------------------------------------------
        self.view_group = QGroupBox("View Settings")

        self.view_mode_label = QLabel("View Mode:")
        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItems(["2D Map Only", "Projected vs PA Only", "Side-by-Side"])

        # Channel selection for the PA plot
        self.pa_chan_label = QLabel("Spectral Channel:")
        self.pa_chan_spin = QSpinBox()
        self.pa_chan_spin.setMinimum(0)
        self.pa_chan_spin.setValue(4)
        self.pa_chan_spin.setToolTip("Index of the wavelength channel for the PA plot")
        self.pa_chan_spin.valueChanged.connect(self._update_wave_display)

        self.wave_val_label = QLabel("λ = N/A")
        self.wave_val_label.setStyleSheet("color: #555; font-style: italic;")

        # --- Plot Options Group -------------------------------------------
        self.options_group = QGroupBox("2D Map Options")

        # 1. Toggles
        self.chk_plot_orbit = QCheckBox("Show Orbit")
        self.chk_plot_orbit.setChecked(True)

        self.chk_colorbar = QCheckBox("Show Colorbar")
        self.chk_colorbar.setChecked(True)

        # 2. Colormap
        self.cmap_label = QLabel("Colormap:")
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(["seismic", "coolwarm", "bwr", "viridis", "plasma", "inferno"])

        # 3. Channel Range (For 2D Map)
        self.chk_full_range = QCheckBox("Use Full Range (Map)")
        self.chk_full_range.setChecked(True)
        self.chk_full_range.toggled.connect(self._update_range_controls)

        self.start_idx_spin = QSpinBox()
        self.start_idx_spin.setPrefix("Start: ")
        self.start_idx_spin.setMinimum(0)

        self.end_idx_spin = QSpinBox()
        self.end_idx_spin.setPrefix("End: ")
        self.end_idx_spin.setMinimum(0)

        # 4. Velocity Shift Slider (-200 to 200 km/s)
        self.shift_label = QLabel("Velocity Shift: 0 km/s")
        self.shift_label.setAlignment(Qt.AlignCenter)
        self.shift_label.setToolTip("Double-click the slider to reset to 0.")

        self.shift_slider = DoubleClickResetSlider(Qt.Horizontal)
        self.shift_slider.setRange(-200, 200)
        self.shift_slider.setValue(0)
        self.shift_slider.setTickPosition(QSlider.TicksBelow)
        self.shift_slider.setTickInterval(50)
        self.shift_slider.valueChanged.connect(self._update_shift_label)

        # --- FOV Group (Symmetric) ----------------------------------------
        self.fov_group = QGroupBox("Field of View (2D Map)")

        self.chk_manual_fov = QCheckBox("Manual Limits")
        self.chk_manual_fov.setChecked(False)
        self.chk_manual_fov.toggled.connect(self._on_manual_fov_toggled)

        # X Limit Slider (0.1 mas to 10 mas)
        self.xlim_slider = QSlider(Qt.Horizontal)
        self.xlim_slider.setRange(10, 1000) # 0.1 to 10.0 mas
        self.xlim_slider.setValue(200)      # Default 2.0 mas
        self.xlim_label = QLabel("H-Zoom: ±2.00 mas")
        self.xlim_slider.valueChanged.connect(self._on_fov_slider_changed)

        # Y Limit Slider
        self.ylim_slider = QSlider(Qt.Horizontal)
        self.ylim_slider.setRange(10, 1000)
        self.ylim_slider.setValue(200)
        self.ylim_label = QLabel("V-Zoom: ±2.00 mas")
        self.ylim_slider.valueChanged.connect(self._on_fov_slider_changed)

        # --- Matplotlib Canvas --------------------------------------------
        self.figure = Figure(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        # --- Logic Connections --------------------------------------------
        self.epoch_combo.currentTextChanged.connect(self.update_baselines)
        self.epoch_combo.currentTextChanged.connect(self.update_files)

        # Initialize
        self.update_baselines()
        self.update_files()
        self._update_range_controls(True)
        self._toggle_fov_sliders(False) # Start disabled
        self._update_wave_display() # Initial update

    def _create_layouts(self):
        # --- Left Side Layout ---
        left_layout = QVBoxLayout()

        # 1. Data Selection
        sel_layout = QFormLayout()
        sel_layout.addRow(self.epoch_label, self.epoch_combo)
        sel_layout.addRow(self.file_label, self.file_combo)
        sel_layout.addRow(self.baseline_label, self.baseline_combo)

        data_vbox = QVBoxLayout()
        data_vbox.addLayout(sel_layout)
        data_vbox.addWidget(self.update_button)
        self.data_group.setLayout(data_vbox)
        left_layout.addWidget(self.data_group)

        # 2. View Settings
        view_layout = QFormLayout()
        view_layout.addRow(self.view_mode_label, self.view_mode_combo)

        # Combine spin and label in one row
        pa_chan_hbox = QHBoxLayout()
        pa_chan_hbox.addWidget(self.pa_chan_spin)
        pa_chan_hbox.addWidget(self.wave_val_label)
        view_layout.addRow(self.pa_chan_label, pa_chan_hbox)

        self.view_group.setLayout(view_layout)
        left_layout.addWidget(self.view_group)

        # 3. Options
        opt_layout = QVBoxLayout()

        # Row 1: Checkboxes
        chk_row = QHBoxLayout()
        chk_row.addWidget(self.chk_plot_orbit)
        chk_row.addWidget(self.chk_colorbar)
        opt_layout.addLayout(chk_row)

        # Row 2: Colormap
        cm_row = QHBoxLayout()
        cm_row.addWidget(self.cmap_label)
        cm_row.addWidget(self.cmap_combo)
        opt_layout.addLayout(cm_row)

        # Row 3: Range
        rng_row = QHBoxLayout()
        rng_row.addWidget(self.chk_full_range)
        rng_row.addWidget(self.start_idx_spin)
        rng_row.addWidget(self.end_idx_spin)
        opt_layout.addLayout(rng_row)

        # Separator line
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        opt_layout.addWidget(line)

        # Row 4: Velocity Shift
        opt_layout.addWidget(self.shift_label)
        opt_layout.addWidget(self.shift_slider)

        self.options_group.setLayout(opt_layout)
        left_layout.addWidget(self.options_group)

        # 4. FOV
        fov_layout = QVBoxLayout()
        fov_layout.addWidget(self.chk_manual_fov)

        fov_layout.addWidget(self.xlim_label)
        fov_layout.addWidget(self.xlim_slider)

        fov_layout.addWidget(self.ylim_label)
        fov_layout.addWidget(self.ylim_slider)

        self.fov_group.setLayout(fov_layout)
        left_layout.addWidget(self.fov_group)

        left_layout.addStretch(1)

        left_widget = QWidget()
        left_widget.setLayout(left_layout)
        left_widget.setMaximumWidth(350) # Keep controls compact

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
    #                            LOGIC HELPERS                           #
    # ------------------------------------------------------------------ #

    def update_baselines(self):
        year = self.epoch_combo.currentText()
        if year:
            opts = [ms.baseline_name(year=year, baseline=i) for i in range(6)]
            self.baseline_combo.clear()
            self.baseline_combo.addItems(opts)

    def update_files(self):
        year = self.epoch_combo.currentText()
        if year in ms.file_num:
            n = ms.file_num[year] + 1
            self.file_combo.clear()
            self.file_combo.addItems([str(i) for i in range(n)])
            self._update_wave_display()

    def _update_shift_label(self):
        val = self.shift_slider.value()
        self.shift_label.setText(f"Velocity Shift: {val} km/s")

    def _update_range_controls(self, full_range):
        self.start_idx_spin.setEnabled(not full_range)
        self.end_idx_spin.setEnabled(not full_range)
        # If toggling on, update the values if possible
        if full_range and self._len_v is not None:
             self.start_idx_spin.setValue(0)
             self.end_idx_spin.setValue(max(0, self._len_v - 1))

    def _update_wave_display(self):
        """Updates the wavelength label based on current file and channel index."""
        year = self.epoch_combo.currentText()
        try:
            filenum = int(self.file_combo.currentText())
            idx = self.pa_chan_spin.value()

            # Get masked wavelength array
            mask = ms.mask_brg(year, filenum)
            waves = ms.WAVE(year, filenum)[mask]

            if 0 <= idx < len(waves):
                w_microns = waves[idx] * 1e6
                self.wave_val_label.setText(f"λ = {w_microns:.5f} μm")
            else:
                self.wave_val_label.setText("λ = out of range")
        except:
            self.wave_val_label.setText("λ = N/A")

    def _on_manual_fov_toggled(self, checked):
        self._toggle_fov_sliders(checked)
        if self._have_plot and checked:
            self.update_plot()

    def _toggle_fov_sliders(self, enabled):
        self.xlim_slider.setEnabled(enabled)
        self.ylim_slider.setEnabled(enabled)
        self.xlim_label.setEnabled(enabled)
        self.ylim_label.setEnabled(enabled)

    def _on_fov_slider_changed(self):
        x_val = self.xlim_slider.value() / self._slider_scale
        y_val = self.ylim_slider.value() / self._slider_scale

        self.xlim_label.setText(f"H-Zoom: ±{x_val:.2f} mas")
        self.ylim_label.setText(f"V-Zoom: ±{y_val:.2f} mas")

        if self._have_plot and self.chk_manual_fov.isChecked():
            if len(self.figure.axes) > 0:
                 ax = self.figure.axes[0]
                 ax.set_xlim(x_val, -x_val)
                 ax.set_ylim(-y_val, y_val)
                 self.canvas.draw_idle()

    def _calculate_orbit_path(self):
        # Parameters
        sma, ecc, inc = 31.091, 0.990, 1.457
        omega = (1.910 + np.pi) % (2*np.pi)
        Omega = (4.603 + np.pi) % (2*np.pi)
        tau, mtot, plx = 0.9894, 14.817, 1.866

        sim = rebound.Simulation()
        sim.units = ('yr', 'AU', 'Msun')

        m1 = 10.29
        sim.add(m=m1)
        sim.add(
            m=mtot-m1, a=sma, e=ecc, inc=inc,
            omega=omega, Omega=Omega,
            T=(1 - tau) * np.sqrt(sma**3/mtot)
        )
        sim.move_to_com()

        t = np.linspace(-10, sim.particles[1].P, 10000)
        x_au = np.empty_like(t)
        y_au = np.empty_like(t)

        for i, ti in enumerate(t):
            sim.integrate(ti)
            p0, p1 = sim.particles[0], sim.particles[1]
            x_au[i] = p1.x - p0.x
            y_au[i] = p1.y - p0.y

        return y_au * plx, x_au * plx

    # ------------------------------------------------------------------ #
    #                              PLOTTING                              #
    # ------------------------------------------------------------------ #
    def update_plot(self):
        year = self.epoch_combo.currentText()
        try:
            filenum = int(self.file_combo.currentText())
        except: return

        # Reset Figure
        self.figure.clear()
        self._cbar = None
        self._have_plot = False

        # Ensure we have valid velocity data first to set limits
        try:
             # Just to get the length
             v_test = ms.pureLinePhotoShift(
                 year=year, filenum=filenum, wave_number=0,
                 return_wave=True, wave=ms.WAVE(year, filenum), showHelp=False
             )
             n_chan = len(v_test)
             self._len_v = n_chan
        except:
             n_chan = 0

        # Update spinbox max limits
        self.start_idx_spin.setMaximum(max(0, n_chan - 1))
        self.end_idx_spin.setMaximum(max(0, n_chan - 1))
        self.pa_chan_spin.setMaximum(max(0, n_chan - 1))

        # If Full Range checked, force values to min/max
        if self.chk_full_range.isChecked():
            self.start_idx_spin.setValue(0)
            self.end_idx_spin.setValue(max(0, n_chan - 1))

        # Update wavelength display just in case
        self._update_wave_display()

        view_mode = self.view_mode_combo.currentText()

        # Layout logic
        if view_mode == "2D Map Only":
            ax = self.figure.add_subplot(111)
            self._plot_2d_map(ax, year, filenum)
        elif view_mode == "Projected vs PA Only":
            ax = self.figure.add_subplot(111)
            self._plot_pa_shift(ax, year, filenum)
        elif view_mode == "Side-by-Side":
            ax1 = self.figure.add_subplot(121)
            ax2 = self.figure.add_subplot(122)
            self._plot_2d_map(ax1, year, filenum)
            self._plot_pa_shift(ax2, year, filenum)
            self.figure.tight_layout()

        self.canvas.draw()
        self._have_plot = True

    # --- PLOT 1: 2D MAP (Original) ------------------------------------
    def _plot_2d_map(self, ax, year, filenum):
        baseline_idx = self.baseline_combo.currentIndex()
        lam_data = ms.WAVE(year, filenum)
        shift_kms = self.shift_slider.value()
        lam = lam_data * (1 - shift_kms / ms.c_km_s)

        ax.set_title(f"Pure-line Shift | {year} | File {filenum}")
        ax.set_xlabel(r"$\Delta$ RA [mas]")
        ax.set_ylabel(r"$\Delta$ Dec [mas]")

        # Continuum photocenter
        try:
            pcx, pcy = ms.orbit(returnPhotocenter=year, plot=False, disableHelp=True)
        except:
            pcx, pcy = 0.0, 0.0

        # Origin & Orbit
        ax.scatter([0], [0], marker="+", c="black", zorder=1000, label="Primary")
        ax.axvline(0, c="k", alpha=0.3, ls="--")
        ax.axhline(0, c="k", alpha=0.3, ls="--")

        if self.chk_plot_orbit.isChecked():
            try:
                ox, oy = self._calculate_orbit_path()
                ax.plot(ox, oy, lw=1.5, alpha=0.6, color="gray", label="Orbit")
            except: pass

        # Velocity coloring
        try:
            v = ms.pureLinePhotoShift(
                year=year, filenum=filenum, wave_number=0,
                return_wave=True, wave=lam, showHelp=False
            )
        except Exception:
            return

        v = np.asarray(v)
        n_chan = len(v)

        # Range Logic
        if self.chk_full_range.isChecked():
            indices = range(n_chan)
        else:
            s = max(0, self.start_idx_spin.value())
            e = min(n_chan-1, self.end_idx_spin.value())
            indices = range(s, e+1)

        # Color Normalization
        valid_v = v[np.isfinite(v)]
        if len(valid_v) > 0:
            vmin_data, vmax_data = valid_v.min(), valid_v.max()
        else:
            vmin_data, vmax_data = -1, 1

        if vmin_data > 0: vmin = -vmax_data
        elif vmax_data < 0: vmax = -vmin_data
        else: vmin, vmax = vmin_data, vmax_data

        cmap_name = self.cmap_combo.currentText()
        cmap = plt.get_cmap(cmap_name)

        try:
            norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        except:
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        colors = [mcolors.to_hex(cmap(norm(val))) if np.isfinite(val) else "#aaa" for val in v]

        # Points
        xs, ys = [], []
        for i in indices:
            if i >= n_chan: continue
            try:
                px, pxe, py, pye = ms.pureLinePhotoShift(
                    year=year, filenum=filenum, baseline=baseline_idx,
                    doFit=True, wave_number=i, plot=True, returnPC=True,
                    wave=lam, showHelp=False
                )
                x, y = px + pcx, py + pcy
                xs.append(x)
                ys.append(y)
                ax.errorbar(
                    x, y, xerr=pxe, yerr=pye, fmt="s", capsize=2,
                    color=colors[i], markeredgecolor="black", alpha=0.9
                )
            except: continue

        # Colorbar
        if self.chk_colorbar.isChecked():
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = self.figure.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("Velocity [km/s]")

        # Limits
        if self.chk_manual_fov.isChecked():
             x_lim = self.xlim_slider.value() / self._slider_scale
             y_lim = self.ylim_slider.value() / self._slider_scale
             ax.set_xlim(x_lim, -x_lim)
             ax.set_ylim(-y_lim, y_lim)
        else:
            if xs:
                all_x = np.array(xs + [0, pcx])
                all_y = np.array(ys + [0, pcy])
                max_r = np.max(np.sqrt(all_x**2 + all_y**2))
                limit = max_r * 1.3 or 1.0
            else:
                limit = 1.0
            ax.set_xlim(limit, -limit)
            ax.set_ylim(-limit, limit)

            # Sync sliders
            self.xlim_slider.blockSignals(True)
            self.ylim_slider.blockSignals(True)
            self.xlim_slider.setValue(int(limit * self._slider_scale))
            self.ylim_slider.setValue(int(limit * self._slider_scale))
            self.xlim_slider.blockSignals(False)
            self.ylim_slider.blockSignals(False)
            self._on_fov_slider_changed()

    # --- PLOT 2: PROJECTED SHIFT vs PA --------------------------------
    def _plot_pa_shift(self, ax, year, filenum):
        wave_number = self.pa_chan_spin.value()

        # 1. Prepare Data Containers
        store_projectedP = []
        store_projectedP_err = []
        store_projectedPtheo = []

        # 2. Extract Basic Data
        try:
            mask = ms.mask_brg(year, filenum)
            # Wavelength for specific channel
            Lam = ms.WAVE(year, filenum)[mask][wave_number]

            UCOORD = ms.Extract(year, filenum, dataType="UCOORD")
            VCOORD = ms.Extract(year, filenum, dataType="VCOORD")
            PA = np.arctan2(UCOORD, VCOORD) # Radians

            # Loop 6 baselines
            for baseline in range(6):
                # Corrected pure line phase
                philine_corr = ms.PHILINE_CORR(
                    year=year, filenum=filenum, baseline=baseline,
                    k_values=None, PhotosphericCorr=False,
                    alpha1_whichSpectrum=None, alpha2_whichSpectrum=None,
                    beta1=None, beta2=None, plot=None, showHelp=False
                )

                # Error
                philine_corr_err = ms.PHILINE(
                    year=year, filenum=filenum, baseline=baseline,
                    PhotosphericCorr=False, returnError=True, showHelp=False
                )

                B_i = np.sqrt(UCOORD**2 + VCOORD**2)[baseline]

                # Projection calculation [mas]
                projectedP = -philine_corr[mask][wave_number] / (2*np.pi) * Lam/B_i * ms.rad2mas
                projectedP_err = np.abs(philine_corr_err[mask][wave_number] / (2*np.pi) * Lam/B_i * ms.rad2mas)

                store_projectedP.append(projectedP)
                store_projectedP_err.append(projectedP_err)

                # Theoretical calculation
                projectedPtheo = -ms.PHI_THEO(
                    year, filenum, baseline, unit="rad"
                )[wave_number] / (2*np.pi) * Lam/B_i * ms.rad2mas
                store_projectedPtheo.append(projectedPtheo)

        except Exception as e:
            ax.text(0.5, 0.5, f"Data Error:\n{str(e)}", ha='center')
            return

        store_projectedP = np.array(store_projectedP)
        store_projectedP_err = np.abs(np.array(store_projectedP_err))
        store_projectedPtheo = np.array(store_projectedPtheo)

        # 3. Fitting
        def cosine_func(x, p, phi):
            return p * np.cos(x - phi)

        # Data for fit (filter NaNs)
        x_data = PA
        y_data = store_projectedP
        y_data_err = store_projectedP_err

        # Colors/Labels
        baseline_labels = np.array([ms.baseline_name(year, b) for b in range(6)])
        baseline_colors = np.array([ms.baseline_colors(b) for b in range(6)])

        nan_mask = np.isfinite(x_data) & np.isfinite(y_data)

        # Fitting Data
        try:
            initial_guess = [3, np.deg2rad(50)]
            bounds = ([0.0, -2*np.pi], [8.0, 2*np.pi])

            popt, pcov = curve_fit(
                cosine_func, x_data[nan_mask], y_data[nan_mask],
                sigma=y_data_err[nan_mask], absolute_sigma=False,
                p0=initial_guess, bounds=bounds, maxfev=50000
            )
            p_fitted, phi_fitted = popt
            perr = np.sqrt(np.diag(pcov))
            p_err, phi_err = perr

            # Fitting Theory
            popt_t, _ = curve_fit(
                cosine_func, PA, store_projectedPtheo,
                p0=initial_guess, bounds=bounds, maxfev=50000
            )
            p_theo, phi_theo = popt_t

            # High res curves
            x_hires = np.linspace(-2*np.pi, 2*np.pi, 200)
            y_hires = cosine_func(x_hires, p_fitted, phi_fitted)
            y_hires_theo = cosine_func(x_hires, p_theo, phi_theo)

            # 4. Plotting
            ax.set_xlabel("Position Angle [deg]")
            ax.set_ylabel("Proj. pure line pc. shift [mas]")
            ax.set_title(f"Shift vs PA | Chan {wave_number}")

            # Plot Theory (Red dashed)
            ax.plot(np.rad2deg(x_hires), y_hires_theo, c='red', ls='--', alpha=0.3, zorder=-1)
            ax.errorbar(np.rad2deg(PA), store_projectedPtheo, fmt='o', c='red', alpha=0.3)
            ax.errorbar(np.rad2deg(PA)+180, -store_projectedPtheo, fmt='o', c='red', alpha=0.3)

            # Plot Fit (Blue)
            ax.plot(np.rad2deg(x_hires), y_hires, c='blue', lw=1.2, label="Fit")

            # Plot Data Points (colored by baseline)
            # Plot both PA and PA+180 (symmetry)
            for i in range(len(baseline_colors)):
                if not nan_mask[i]: continue

                # Original point
                ax.errorbar(
                    np.rad2deg(x_data[i]), y_data[i], yerr=y_data_err[i],
                    fmt='s', color=baseline_colors[i], markeredgecolor='black',
                    label=baseline_labels[i]
                )
                # Symmetric point
                ax.errorbar(
                    np.rad2deg(x_data[i])+180, -y_data[i], yerr=y_data_err[i],
                    fmt='s', color=baseline_colors[i], markeredgecolor='black'
                )

            # Text Box
            param_text = (
                rf"$p = {p_fitted:.2f} \pm {p_err:.2f}\,\mathrm{{mas}}$" "\n"
                rf"$\varphi = {np.rad2deg(phi_fitted):.1f} \pm {np.rad2deg(phi_err):.1f}\,^\circ$"
            )
            ax.text(0.02, 0.98, param_text, transform=ax.transAxes, ha="left", va="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

            # Limits specific to year (from snippet)
            if year == "2019_red":
                ax.set_ylim(-1, 1)
            else:
                ax.set_ylim(-6.5, 5)

            # Handle legend - remove duplicates
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), ncol=2, loc="lower left", fontsize="small")

        except Exception as e:
            ax.text(0.5, 0.5, f"Fit Error:\n{str(e)}", ha='center')
