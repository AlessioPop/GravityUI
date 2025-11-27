# GravityUI_main.py

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget

from visibility_tab import VisibilityTab
from normalized_flux_tab import NormalizedFluxTab
from differential_phase_tab import DifferentialPhaseTab
from photocenter_shift_tab import PhotocenterShiftTab
from uv_coverage_tab import UVCoverageTab  # Import new tab
from style_sheet import apply_style


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("GravityUI")
        self.resize(1200, 750)

        # Central tab widget
        self.tabs = QTabWidget()
        self.tabs.setObjectName("MainTabs")
        self.setCentralWidget(self.tabs)

        # --- Add Tabs ---
        
        # 1. U-V Coverage (New)
        self.uv_tab = UVCoverageTab()
        self.tabs.addTab(self.uv_tab, "U-V Coverage")

        # 2. Visibility
        self.vis_tab = VisibilityTab()
        self.tabs.addTab(self.vis_tab, "Visibility")

        # 3. Differential Phase
        self.diffphi_tab = DifferentialPhaseTab()
        self.tabs.addTab(self.diffphi_tab, "Diff. Phase")

        # 4. Normalized Flux
        self.normflux_tab = NormalizedFluxTab()
        self.tabs.addTab(self.normflux_tab, "Norm. Flux")

        # 5. Photocenter Shift
        self.photocenter_tab = PhotocenterShiftTab()
        self.tabs.addTab(self.photocenter_tab, "Photocenter Shift")

def main():
    app = QApplication(sys.argv)

    # Apply global style
    apply_style(app)

    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()