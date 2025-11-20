# style_sheet.py
#
# Central QSS + palette for GravityUI.
# Minimal, scientific, light theme with card-style panels.

from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtCore import Qt


MAIN_QSS = """
/* -------- General window background -------- */
QMainWindow {
    background-color: #f3f5fa;
}

/* -------- Tab widget / tab bar -------- */
QTabWidget::pane {
    border: none;
    background: transparent;
}

QTabBar::tab {
    background: #e5e9f2;
    color: #5c6475;
    border-radius: 14px;
    padding: 6px 16px;
    margin-right: 4px;
    font-weight: 500;
}

QTabBar::tab:selected {
    background: #ffffff;
    color: #111827;
}

QTabBar::tab:hover:!selected {
    background: #edf0f7;
}

/* -------- Left side panel "card" -------- */
/* You will set objectName("SidePanel") on the left column container. */
QWidget#SidePanel {
    background-color: #ffffff;
    border-radius: 18px;
    padding: 16px;
}

/* Optional: a white card for parameter groups, if you set objectName("ParamCard") */
QGroupBox#ParamCard {
    background-color: #ffffff;
    border-radius: 14px;
    padding: 10px 12px 12px 12px;
    border: 1px solid #e2e5f0;
    margin-top: 8px;
}

/* -------- Group boxes (headings on the left) -------- */
QGroupBox {
    border: none;
    margin-top: 12px;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 0;
    padding: 0 0 4px 0;
    font-weight: 600;
    color: #111827;
}

/* -------- Labels -------- */
QLabel {
    color: #111827;
}

/* Subtitle / helper text (if you want to use objectName) */
QLabel#MutedLabel {
    color: #8a92a6;
    font-size: 9pt;
}

/* -------- Buttons -------- */
QPushButton {
    background-color: #f3f4ff;
    color: #111827;
    border-radius: 12px;
    border: 1px solid #d5d9e6;
    padding: 6px 14px;
}

QPushButton:hover {
    background-color: #e9ecff;
}

QPushButton:pressed {
    background-color: #d7dcff;
}

QPushButton:disabled {
    background-color: #eef0f6;
    color: #b7bcc9;
    border-color: #e0e3ed;
}

/* -------- Combo boxes -------- */
QComboBox {
    background-color: #ffffff;
    border-radius: 10px;
    border: 1px solid #d5d9e6;
    padding: 4px 10px;
}

QComboBox::drop-down {
    border: none;
    width: 20px;
}

QComboBox::down-arrow {
    image: none;  /* use native arrow */
}

/* -------- Checkboxes -------- */
QCheckBox {
    spacing: 6px;
    color: #111827;
}

QCheckBox::indicator {
    width: 14px;
    height: 14px;
    border-radius: 4px;
    border: 1px solid #c6cad8;
    background: #ffffff;
}

QCheckBox::indicator:hover {
    border-color: #ffb468;
}

QCheckBox::indicator:checked {
    background-color: #ff9b42;
    border-color: #ff9b42;
}

/* -------- Scroll areas (if any) -------- */
QScrollArea {
    border: none;
    background: transparent;
}

/* -------- Matplotlib toolbar -------- */
QToolBar {
    background-color: transparent;
    border: none;
}

/* Canvas is a QWidget; give it a clean white card look if desired. */
QWidget#PlotCanvas {
    background-color: #ffffff;
    border-radius: 18px;
}
"""


def apply_style(app):
    """Apply Fusion style, a light palette, and the QSS."""
    app.setStyle("Fusion")

    palette = QPalette()
    palette.setColor(QPalette.Window, QColor("#f3f5fa"))
    palette.setColor(QPalette.Base, QColor("#ffffff"))
    palette.setColor(QPalette.AlternateBase, QColor("#f3f5fa"))
    palette.setColor(QPalette.Text, QColor("#111827"))
    palette.setColor(QPalette.Button, QColor("#f3f4ff"))
    palette.setColor(QPalette.ButtonText, QColor("#111827"))
    palette.setColor(QPalette.Highlight, QColor("#ff9b42"))
    palette.setColor(QPalette.HighlightedText, QColor("#ffffff"))
    palette.setColor(QPalette.Link, QColor("#ff9b42"))

    app.setPalette(palette)
    app.setStyleSheet(MAIN_QSS)
