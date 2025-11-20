# GravityUI_main.py

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget

from visibility_tab import VisibilityTab
from normalized_flux_tab import NormalizedFluxTab
from differential_phase_tab import DifferentialPhaseTab
from style_sheet import apply_style


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("GravityUI")
        self.resize(1200, 700)

        # Central tab widget
        self.tabs = QTabWidget()
        self.tabs.setObjectName("MainTabs")
        self.setCentralWidget(self.tabs)

        # Tabs
        self.vis_tab = VisibilityTab()
        self.tabs.addTab(self.vis_tab, "Visibility Amplitude")

        self.diffphi_tab = DifferentialPhaseTab()
        self.tabs.addTab(self.diffphi_tab, "Differential Phase")

        self.normflux_tab = NormalizedFluxTab()
        self.tabs.addTab(self.normflux_tab, "Normalized Flux")


def main():
    app = QApplication(sys.argv)

    # Apply global style (Fusion+QSS)
    apply_style(app)

    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
