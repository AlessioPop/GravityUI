# style_sheet.py
#
# Central QSS + palette for GravityUI.
# Theme: "Classic Professional" - Clean, light, standard UI with high readability.

from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtCore import Qt

# --- Color Palette (Light/Professional) ---
C_BG_MAIN       = "#F5F5F5"  # Main window background (Off-White)
C_BG_PANEL      = "#FFFFFF"  # Panels/Groupboxes (White)
C_BG_INPUT      = "#FFFFFF"  # Input fields
C_BG_BUTTON     = "#E0E0E0"  # Button gradient start
C_BG_HOVER      = "#D5D5D5"  # Hover state

# Text
C_TEXT_PRI      = "#000000"  # Black
C_TEXT_SEC      = "#555555"  # Dark Grey
C_TEXT_DIS      = "#A0A0A0"  # Disabled Grey

# Accents
C_ACCENT        = "#0078D7"  # Standard Windows Blue
C_ACCENT_HOVER  = "#005A9E"  # Darker Blue
C_BORDER        = "#A0A0A0"  # Standard borders

MAIN_QSS = f"""
/* ==========================================================================
   GLOBAL RESET & MAIN WINDOW
   ========================================================================== */
QMainWindow {{
    background-color: {C_BG_MAIN};
    color: {C_TEXT_PRI};
}}

QWidget {{
    font-family: "Segoe UI", "Roboto", "Helvetica", sans-serif;
    font-size: 10pt;
    color: {C_TEXT_PRI};
}}

/* ==========================================================================
   TAB WIDGET
   ========================================================================== */
QTabWidget::pane {{
    border: 1px solid {C_BORDER};
    background-color: {C_BG_PANEL};
    border-radius: 3px;
    top: -1px;
}}

QTabBar::tab {{
    background: {C_BG_BUTTON};
    color: {C_TEXT_PRI};
    border: 1px solid {C_BORDER};
    border-bottom: 1px solid {C_BORDER};
    border-top-left-radius: 3px;
    border-top-right-radius: 3px;
    padding: 6px 12px;
    margin-right: 2px;
}}

QTabBar::tab:selected {{
    background: {C_BG_PANEL};
    border-bottom: 1px solid {C_BG_PANEL}; /* Blend into pane */
    font-weight: bold;
}}

QTabBar::tab:hover:!selected {{
    background: {C_BG_HOVER};
}}

/* ==========================================================================
   CONTAINERS
   ========================================================================== */
QGroupBox {{
    background-color: {C_BG_PANEL};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    margin-top: 20px;
    padding-top: 16px;
    font-weight: bold;
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 2px 6px;
    left: 8px;
    color: {C_TEXT_PRI};
}}

/* ==========================================================================
   INPUTS (Standard Look)
   ========================================================================== */
QComboBox, QSpinBox, QDoubleSpinBox {{
    background-color: {C_BG_INPUT};
    border: 1px solid {C_BORDER};
    border-radius: 3px;
    padding: 4px 8px;
    color: {C_TEXT_PRI};
    min-height: 20px;
}}

QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover {{
    border: 1px solid {C_ACCENT};
}}

/* We let Qt draw the default arrow, just adding padding so it fits */
QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 24px;
    border-left-width: 1px;
    border-left-color: {C_BORDER};
    border-left-style: solid;
    border-top-right-radius: 3px;
    border-bottom-right-radius: 3px;
    background: {C_BG_BUTTON};
}}

QComboBox::down-arrow {{
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid black; /* Simple black triangle */
    margin-top: 2px;
    margin-right: 2px;
}}

/* Spinbox Buttons */
QSpinBox::up-button, QDoubleSpinBox::up-button {{
    subcontrol-origin: border;
    subcontrol-position: top right;
    width: 18px;
    background: {C_BG_BUTTON};
    border-left: 1px solid {C_BORDER};
    border-bottom: 1px solid {C_BORDER};
}}

QSpinBox::down-button, QDoubleSpinBox::down-button {{
    subcontrol-origin: border;
    subcontrol-position: bottom right;
    width: 18px;
    background: {C_BG_BUTTON};
    border-left: 1px solid {C_BORDER};
    border-top: 0px solid {C_BORDER};
}}

/* Simple black arrows for spinbox */
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
    width: 0; height: 0;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-bottom: 5px solid black;
}}

QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
    width: 0; height: 0;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid black;
}}

/* ==========================================================================
   BUTTONS
   ========================================================================== */
QPushButton {{
    background-color: {C_BG_BUTTON};
    color: {C_TEXT_PRI};
    border: 1px solid {C_BORDER};
    border-radius: 3px;
    padding: 6px 16px;
}}

QPushButton:hover {{
    background-color: {C_BG_HOVER};
    border: 1px solid {C_ACCENT};
}}

QPushButton:pressed {{
    background-color: {C_ACCENT};
    color: white;
}}

/* ==========================================================================
   CHECKBOX
   ========================================================================== */
QCheckBox {{
    spacing: 6px;
}}

QCheckBox::indicator {{
    width: 14px;
    height: 14px;
    border: 1px solid {C_BORDER};
    background: {C_BG_INPUT};
    border-radius: 2px;
}}

QCheckBox::indicator:hover {{
    border-color: {C_ACCENT};
}}

QCheckBox::indicator:checked {{
    background-color: {C_ACCENT};
    border-color: {C_ACCENT};
    image: none; /* Using solid color for simplicity */
}}

/* ==========================================================================
   SLIDERS
   ========================================================================== */
QSlider::groove:horizontal {{
    border: 1px solid {C_BORDER};
    height: 6px;
    background: {C_BG_BUTTON};
    margin: 2px 0;
    border-radius: 3px;
}}

QSlider::handle:horizontal {{
    background: {C_ACCENT};
    border: 1px solid {C_ACCENT};
    width: 12px;
    height: 12px;
    margin: -4px 0;
    border-radius: 6px;
}}

QSlider::handle:horizontal:hover {{
    background: {C_ACCENT_HOVER};
}}

/* ==========================================================================
   MISC
   ========================================================================== */
QScrollArea {{
    border: none;
    background: transparent;
}}

QLabel {{
    color: {C_TEXT_PRI};
}}

/* Matplotlib Canvas Frame */
QWidget > FigureCanvasQTAgg {{
    border: 1px solid {C_BORDER};
    background-color: {C_BG_PANEL};
}}
"""

def apply_style(app):
    """Apply the Classic Professional (Light) theme."""
    app.setStyle("Fusion")

    palette = QPalette()
    palette.setColor(QPalette.Window,          QColor(C_BG_MAIN))
    palette.setColor(QPalette.WindowText,      QColor(C_TEXT_PRI))
    palette.setColor(QPalette.Base,            QColor(C_BG_INPUT))
    palette.setColor(QPalette.AlternateBase,   QColor(C_BG_PANEL))
    palette.setColor(QPalette.Text,            QColor(C_TEXT_PRI))
    palette.setColor(QPalette.Button,          QColor(C_BG_BUTTON))
    palette.setColor(QPalette.ButtonText,      QColor(C_TEXT_PRI))
    palette.setColor(QPalette.Highlight,       QColor(C_ACCENT))
    palette.setColor(QPalette.HighlightedText, QColor("#FFFFFF"))

    # Disabled colors
    palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(C_TEXT_DIS))
    palette.setColor(QPalette.Disabled, QPalette.Text,       QColor(C_TEXT_DIS))
    palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(C_TEXT_DIS))

    app.setPalette(palette)
    app.setStyleSheet(MAIN_QSS)
