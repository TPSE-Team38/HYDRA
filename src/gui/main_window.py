import os
from datetime import datetime
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QLabel, QPushButton, QFileDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox,
    QLineEdit, QSplitter, QTextEdit, QMessageBox, QScrollArea, QApplication, QComboBox, QTextBrowser
)
from PySide6.QtCore import Qt, QThreadPool
from PySide6.QtGui import QIcon
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from .PlotFullscreenWindow import PlotFullscreenWindow
from .plot_widget import PlotWidget
from .controller import AnalysisController
from .protein_row import ProteinInputRow
from src.models import AnalysisConfig
from src.Calculations import tau, peclet, get_z_vals
import src.gui.accessibility_colors as accessibility_colors
from PySide6.QtWidgets import QLabel, QToolTip
from PySide6.QtGui import QCursor,QBitmap

from ..parallization import LoadMS1Worker


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.plot_toolbar = None
        self.accessibility_win = None
        self.plot_window = None

        self.setWindowTitle("HYDRA")
        self.resize(1200, 800)
        self.reset_btn,self.abort_remasking_btn,self.continue_remasking_btn =None,None,None
        self.icon=QIcon("./assets/HYDRA-logo.png")
        self.setWindowIcon(self.icon)
        self.ms1_path = None
        self.protein_rows = []
        self.analysis_results = []  # list of results
        self.current_result_index = 0  # which protein is shown

        self.plot = PlotWidget()

        self.controller = AnalysisController(self.plot, self)
        self._build_ui()
        self.info_box.setOpenLinks(False)  # Prevents clicking from doing anything
        self.info_box.highlighted.connect(self.show_definition_tooltip)
    # ================= UI =================
    def _build_ui(self):
        central = QWidget()
        main_layout = QVBoxLayout(central)

        # ---------- TOP BAR ----------
        top_bar = QHBoxLayout()

        self.load_btn = QPushButton("Upload ms1 file")
        self.load_btn.clicked.connect(self.load_file)

        self.file_status = QLabel("No ms1 file loaded")
        self.file_status.setStyleSheet("color: red;")

        title = QLabel("HYDRA")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 21px; font-weight: bold;")

        self.accessibility_btn = QPushButton("Accessibility")
        self.accessibility_btn.clicked.connect(self.open_accessibility_settings)
        top_bar.addWidget(self.load_btn)
        top_bar.addWidget(self.file_status)
        top_bar.addStretch()
        top_bar.addWidget(title)
        top_bar.addStretch()
        top_bar.addWidget(self.accessibility_btn)

        main_layout.addLayout(top_bar)

        # ---------- MAIN SPLITTER ----------
        main_splitter = QSplitter(Qt.Orientation.Vertical)
        main_layout.addWidget(main_splitter)

        # ================= UPPER PANEL =================
        upper_widget = QWidget()
        upper_layout = QVBoxLayout(upper_widget)

        upper_splitter = QSplitter(Qt.Orientation.Horizontal)

        # ---------- Protein Input ----------
        protein_box = QGroupBox("Protein Input")
        protein_layout = QVBoxLayout(protein_box)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        scroll_content = QWidget()
        self.protein_list_layout = QVBoxLayout(scroll_content)
        scroll.setWidget(scroll_content)

        protein_layout.addWidget(scroll)

        add_protein_btn = QPushButton("+ Add Protein")
        add_protein_btn.clicked.connect(self.add_protein_row)
        protein_layout.addWidget(add_protein_btn)

        self.add_protein_row()

        # ---------- Parameters ----------
        param_box = QGroupBox("Parameters")
        param_layout = QGridLayout(param_box)

        self.temp_input = QLineEdit("22")
        self.viscosity_input = QLineEdit("0.9544e-3")
        self.radius_input = QLineEdit("125")
        self.length_input = QLineEdit("114.5")
        self.flow_input = QLineEdit("20")



        params = [
            ("Temperature (°C)", self.temp_input),
            ("Viscosity (kg m-1 s-1)", self.viscosity_input),
            ("Capillary radius (µm)", self.radius_input),
            ("Capillary length (cm)", self.length_input),
            ("Flow rate (µL/min)", self.flow_input),
        ]

        for row, (label, widget) in enumerate(params):
            param_layout.addWidget(QLabel(label), row, 0)
            param_layout.addWidget(widget, row, 1)

        # ---------- Tau / Péclet ----------
        tau_box = QGroupBox("Tau and Péclet")
        tau_layout = QVBoxLayout(tau_box)

        self.estimated_rh_input = QLineEdit("2.10")


        self.calc_tau_btn = QPushButton("Calculate Tau and Péclet")
        self.calc_tau_btn.setStyleSheet("font-weight: bold;")
        self.calc_tau_btn.clicked.connect(self.calculate_tau_peclet)

        self.tau_out = QLineEdit()
        self.tau_out.setReadOnly(True)

        self.peclet_out = QLineEdit()
        self.peclet_out.setReadOnly(True)

        tau_layout.addWidget(QLabel("Estimated Rh (nm)"))
        tau_layout.addWidget(self.estimated_rh_input)
        tau_layout.addWidget(self.calc_tau_btn)
        tau_layout.addWidget(QLabel("Tau (τ)"))
        tau_layout.addWidget(self.tau_out)
        tau_layout.addWidget(QLabel("Péclet"))
        tau_layout.addWidget(self.peclet_out)

        upper_splitter.addWidget(protein_box)
        upper_splitter.addWidget(param_box)
        upper_splitter.addWidget(tau_box)

        upper_layout.addWidget(upper_splitter)

        # ---------- Analyse Button ----------
        self.analyse_btn = QPushButton("Analyse Data")
        self.analyse_btn.setStyleSheet("font-size: 16px; font-weight: bold;")
        self.analyse_btn.clicked.connect(self.run)
        upper_layout.addWidget(self.analyse_btn, alignment=Qt.AlignmentFlag.AlignCenter)

        # ---------- Export Button ----------
        self.export_btn = QPushButton("Export Results as PDF")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self.export_pdf)
        upper_layout.addWidget(self.export_btn, alignment=Qt.AlignmentFlag.AlignCenter)

        main_splitter.addWidget(upper_widget)

        # ================= LOWER PANEL =================

        lower_splitter = QSplitter(Qt.Orientation.Horizontal)

        # ---------- Plot ----------
        plot_container = QWidget()
        self.plot_layout = QVBoxLayout(plot_container)

        # Remove margins so the plot fits snugly
        self.plot_layout.setContentsMargins(0, 0, 0, 0)

        self.plot_toolbar = NavigationToolbar(self.plot, self)

        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Results Plot"))

        open_fullscreen_btn = QPushButton("View Plot in Full Screen")
        open_fullscreen_btn.clicked.connect(self.open_plot_fullscreen)

        header_layout.addStretch()
        header_layout.addWidget(open_fullscreen_btn)

        self.plot_layout.addLayout(header_layout)
        self.plot_layout.addWidget(self.plot_toolbar)
        self.plot_layout.addWidget(self.plot)

        lower_splitter.addWidget(plot_container)

        # info box
        self.info_box = QTextBrowser()
        self.info_box.setReadOnly(True)
        self.info_box.setMinimumWidth(250)
        self.info_box.setMaximumWidth(400)
        self.info_box.setPlaceholderText("Analysis results will appear here")

        lower_splitter.addWidget(self.info_box)

        lower_splitter.setStretchFactor(0, 3)
        lower_splitter.setStretchFactor(1, 1)

        main_splitter.addWidget(lower_splitter)

        main_splitter.setStretchFactor(0, 0)
        main_splitter.setStretchFactor(1, 1)

        self.setCentralWidget(central)

        # Navigation buttons

        nav_layout = QHBoxLayout()

        self.prev_btn = QPushButton("< Previous")
        self.abort_remasking_btn = QPushButton("Abort Remasking Ø")
        self.continue_remasking_btn = QPushButton("Confirm Remasking ✓")
        self.next_btn = QPushButton("Next >")

        self.prev_btn.clicked.connect(self.show_previous)
        self.next_btn.clicked.connect(self.show_next)
        self.prev_btn.setEnabled(False)
        self.next_btn.setEnabled(False)
        self.continue_remasking_btn.setEnabled(False)
        self.abort_remasking_btn.setEnabled(False)
        # ---------- Reset Button ----------
        self.reset_btn = QPushButton("Reset Masking")
        self.reset_btn.setEnabled(False)
        nav_layout.addWidget(self.prev_btn)
        # nav_layout.addStretch()
        nav_layout.addWidget(self.abort_remasking_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.reset_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.continue_remasking_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.next_btn)

        self.navigation_layout = nav_layout

        main_layout.addLayout(nav_layout)

    def show_definition_tooltip(self, qurl):
        """Catches the hover over 'Tau' and shows the definition."""
        definition_text = qurl.toString()

        if definition_text:
            # Show the popup at the mouse cursor
            QToolTip.showText(QCursor.pos(), definition_text)
        else:
            # Hide it when moving the mouse away
            QToolTip.hideText()

    # ================= ACTIONS =================
    def open_accessibility_settings(self,checked):
        self.accessibility_win=Accessibility_Window(self,accessibility_colors.current_mode)
        self.accessibility_win.show()

    def restore_plot_to_main(self):
        # Detach from fullscreen
        self.plot.setParent(None)
        self.plot_toolbar.setParent(None)

        self.prev_btn.setParent(None)
        self.next_btn.setParent(None)
        self.reset_btn.setParent(None)
        self.abort_remasking_btn.setParent(None)
        self.continue_remasking_btn.setParent(None)

        # Restore plot
        self.plot_layout.addWidget(self.plot_toolbar)
        self.plot_layout.addWidget(self.plot)

        # Restore buttons to main layout (adjust layout name if different)
        self.navigation_layout.addWidget(self.prev_btn)
        self.navigation_layout.addWidget(self.abort_remasking_btn)
        self.navigation_layout.addWidget(self.reset_btn)
        self.navigation_layout.addWidget(self.continue_remasking_btn)
        self.navigation_layout.addWidget(self.next_btn)

        self.show_current_result()



    def open_plot_fullscreen(self):

        if self.plot.stored_show_eic_args is None:
            QMessageBox.warning(self, "No plot", "No plot available yet.")
            return

        if (
                self.plot.curr_res is not None
                and self.plot.curr_res.currently_remasking
        ):
            QMessageBox.warning(
                self,
                "Remasking in progress",
                "Please confirm or abort remasking before opening fullscreen view."
            )
            return


        # Detach from main window
        self.plot.setParent(None)
        self.plot_toolbar.setParent(None)
        # Detach buttons from main window
        self.prev_btn.setParent(None)
        self.next_btn.setParent(None)
        self.reset_btn.setParent(None)
        self.abort_remasking_btn.setParent(None)
        self.continue_remasking_btn.setParent(None)

        self.plot_window = PlotFullscreenWindow(
            parent=self,
            plot=self.plot,
            toolbar=self.plot_toolbar,
            on_close=self.restore_plot_to_main
        )

        self.plot_window.setWindowTitle("Results Plot")
        self.plot_window.setWindowIcon(self.icon)
        self.plot_window.resize(1400, 900)

        central = QWidget()
        layout = QVBoxLayout(central)

        # Plot + toolbar
        layout.addWidget(self.plot_toolbar)
        layout.addWidget(self.plot)

        # ---- reuse main window buttons ----
        nav_layout = QHBoxLayout()
        nav_layout.addWidget(self.prev_btn)
        nav_layout.addWidget(self.abort_remasking_btn)
        nav_layout.addWidget(self.reset_btn)
        nav_layout.addWidget(self.continue_remasking_btn)
        nav_layout.addWidget(self.next_btn)
        layout.addLayout(nav_layout)

        self.plot_window.setCentralWidget(central)
        self.show_current_result()
        self.plot_window.showMaximized()


    def reset_masking(self):
        self.show_current_result()


    def show_recalculated_fit(self,masked_y, fitted_y, r2, t_R, sigma, D, R_h, t, p):
        result=self.analysis_results[self.current_result_index]

        r2_value = f"{result.r2}"
        if result.r2 >= 0.95:
            r2_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {r2_value} ✔</span>'
        elif 0.95 > result.r2 >= 0.9:
            r2_colored = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {r2_value} ⚠️</span>'
        else:
            r2_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {r2_value} ❌</span>'

        tau_value = f"{result.t}"
        if result.t >= 1.25:
            tau_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {tau_value} ✔</span>'
        elif result.t > 0.37:
            tau_colored = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {tau_value} ⚠️</span>'
        else:
            tau_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {tau_value} ❌</span>'

        peclet_value = f"{result.p}"
        if result.p > 70:
            peclet_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {peclet_value} ✔</span>'
        else:
            peclet_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {peclet_value} ❌</span>'

        r2_value_remasked = f"{r2}"
        if r2 >= 0.95:
            r2_colored_remasked = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {r2_value_remasked} ✔</span>'
        elif 0.95 > r2 >= 0.9:
            r2_colored_remasked = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {r2_value_remasked} ⚠️</span>'
        else:
            r2_colored_remasked = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {r2_value_remasked} ❌</span>'

        tau_value_remasked = f"{t}"
        if t >= 1.25:
            tau_colored_remasked = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {tau_value_remasked} ✔</span>'
        elif t > 0.37:
            tau_colored_remasked = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {tau_value_remasked} ⚠️</span>'
        else:
            tau_colored_remasked = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {tau_value_remasked} ❌</span>'

        peclet_value_remasked = f"{p}"
        if p > 70:
            peclet_colored_remasked = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {peclet_value_remasked} ✔</span>'
        else:
            peclet_colored_remasked = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {peclet_value_remasked} ❌</span>'

        z_vals = get_z_vals(result.charge_state, result.charge_range)
        z_text = ", ".join(str(int(z)) for z in z_vals)

        tau_definition = f"Tau is the time it takes for the protein to reach 95% of its maximum fluorescence intensity"
        peclet_definition = f"Péclet number is the ratio of the protein's diffusion coefficient to the flow rate"
        plain_text_style = "text-decoration: none; color: inherit;"

        self.update_info(
            f"Protein {self.current_result_index + 1} / {len(self.analysis_results)}<br>"
            f"{self.protein_rows[self.current_result_index].proteinName.text()}<br>"
            f"m/z: {result.protein_mz} range {result.mz_window}<br>"
            f"Charge {result.charge_state} with range {result.charge_range} includes charges:<br>{z_text}<br><br>"
            "Previous Fit:<br> <br>"
            f"t_R: {result.tR:.2f} s<br>"
            f"σ: {result.sigma:.3e}<br>"
            f"R²: {r2_colored}<br>"
            f"<b>R_h: {result.Rh:.3e} m </b><br>"
            f"D: {result.D:.3e} m²/s<br>"
            f"<a href='{tau_definition}' style='{plain_text_style}'>Tau</a>: {tau_colored}<br>"
            f"<a href='{peclet_definition}' style='{plain_text_style}'>Péclet</a>: {peclet_colored}<br> <br>"
            "Recalculated Fit:<br>"
            f"t_R: {t_R:.2f} s<br>"
            f"σ: {sigma:.3e} s<br>"
            f"R²: {r2_colored_remasked}<br>"
            f"<b>R_h: {R_h:.3e} m</b><br>"
            f"D: {D:.3e} m²/s<br>"
            f"<a href='{tau_definition}' style='{plain_text_style}'>Tau</a>: {tau_colored_remasked}<br>"
            f"<a href='{peclet_definition}' style='{plain_text_style}'>Péclet</a>: {peclet_colored_remasked}<br>"
        )

    def add_protein_row(self):
        row = ProteinInputRow(parent_window=self)
        self.protein_rows.append(row)
        self.protein_list_layout.addWidget(row)

    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select ms1 file", "", "*.ms1")
        if path:
            self.file_status.setText("Loading...")
            self.file_status.setStyleSheet("color: orange; font-weight: bold;")
            self.load_btn.setEnabled(False)
            QApplication.processEvents()
            worker=LoadMS1Worker(self.controller.load_ms1_once,path)
            worker.signals.finished.connect(self.on_ms1_loaded)
            worker.signals.error.connect(self.on_ms1_error)
            QThreadPool.globalInstance().start(worker)
            # self.threadpool.globalInstance().start(self.controller.load_ms1_once(path))
            # self.controller.load_ms1_once(path)
            # self.ms1_path = path
            # self.intermediate_path=path
            # self.file_status.setText("✔ ms1 file loaded")
            # self.file_status.setStyleSheet("color: green;")

    def on_ms1_loaded(self,path):
        self.ms1_path = path
        self.load_btn.setEnabled(True)
        file_name = os.path.basename(path)
        self.file_status.setText(f"✔ ms1 file loaded: {self.ms1_path}")
        self.file_status.setStyleSheet("color: green;")

    def on_ms1_error(self, msg):
        self.file_status.setText("❌✔ Failed to load ms1")
        self.file_status.setStyleSheet("color: red;")
        self.load_btn.setEnabled(True)
        print(msg)

    def run(self):
        if not self.ms1_path:
            QMessageBox.warning(self, "Error", "Please load an MS1 file first.")
            return

        if not self.protein_rows:
            QMessageBox.warning(self, "Error", "Please add at least one protein.")
            return

        self.analysis_results.clear()

        valid_rows = [row for row in self.protein_rows if self.is_row_valid(row)]

        if not valid_rows:
            QMessageBox.warning(
                self,
                "No valid protein",
                "Please enter at least one complete protein row."
            )
            return
        self.analyse_btn.setText("Analysing...")
        self.analyse_btn.setEnabled(True)
        QApplication.processEvents()
        try:
            for i, row in enumerate(valid_rows, start=1):
                config = AnalysisConfig(
                    ms1_path=self.ms1_path,
                    protein_name=row.proteinName.text(),
                    protein_mz=float(row.mz.text()),
                    mz_window=float(row.range.text()),
                    charge_state=int(row.charge.text()),
                    charge_range=int(row.charge_range.text()),
                    temperature=float(self.temp_input.text()),
                    viscosity=float(self.viscosity_input.text()),
                    capillary_radius=float(self.radius_input.text()),
                    capillary_length=float(self.length_input.text()),
                    flow_rate=float(self.flow_input.text()),
                )

                result = self.controller.run(config,self.reset_btn,self.abort_remasking_btn,self.continue_remasking_btn,self.show_current_result, store=False)
                if result is None:
                    QMessageBox.critical(self, "Analysis failed", "Analysis returned no result.")
                    return
                self.analysis_results.append(result)
            self.analyse_btn.setText("Analyse Data")
            self.analyse_btn.setEnabled(True)
            self.export_btn.setEnabled(True)
            self.current_result_index = 0
            self.show_current_result()

        except ValueError:
            QMessageBox.critical(self, "Invalid input", "Check protein or parameter values.")

    def show_current_result(self):
        if not self.analysis_results:
            return

        result = self.analysis_results[self.current_result_index]


        self.plot.show_eic(result, self.reset_btn,self.show_current_result,self.abort_remasking_btn,self.continue_remasking_btn,self.show_recalculated_fit,config=AnalysisConfig(
            ms1_path=self.ms1_path,
            protein_name=result.protein_name,
            protein_mz=float(result.protein_mz),
            mz_window=float(result.mz_window),
            charge_state=int(result.charge_state),
            charge_range=int(result.charge_range),
            temperature=float(self.temp_input.text()),
            viscosity=float(self.viscosity_input.text()),
            capillary_radius=float(self.radius_input.text()),
            capillary_length=float(self.length_input.text()),
            flow_rate=float(self.flow_input.text())
        ))

        r2_value = f"{result.r2}"
        if result.r2 >= 0.95:
            r2_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {r2_value} ✔</span>'
        elif 0.95 > result.r2 >= 0.9:
            r2_colored = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {r2_value} ⚠️</span>'
        else:
            r2_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {r2_value} ❌</span>'

        tau_value = f"{result.t}"
        if result.t >= 1.25:
            tau_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {tau_value} ✔</span>'
        elif result.t > 0.37:
            tau_colored = f'<span style = "color:{accessibility_colors.YellowBox_Border}"> {tau_value} ⚠️</span>'
        else:
            tau_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {tau_value} ❌</span>'

        peclet_value = f"{result.p}"
        if result.p > 70:
            peclet_colored = f'<span style = "color:{accessibility_colors.GreenBox_Border}"> {peclet_value} ✔</span>'
        else:
            peclet_colored = f'<span style = "color:{accessibility_colors.RedBox_Border}"> {peclet_value} ❌</span>'

        z_vals = get_z_vals(result.charge_state, result.charge_range)
        z_text = ", ".join(str(int(z)) for z in z_vals)

        tau_definition = f"Tau is the time it takes for the protein to reach 95% of its maximum fluorescence intensity"
        peclet_definition = f"Péclet number is the ratio of the protein's diffusion coefficient to the flow rate"
        plain_text_style = "text-decoration: none; color: inherit;"
        self.update_info(

            f"Protein {self.current_result_index + 1} / {len(self.analysis_results)}<br>"
            f"{self.protein_rows[self.current_result_index].proteinName.text()}<br>"
            f"m/z: {result.protein_mz} range {result.mz_window}<br>"
            f"Charge {result.charge_state} with range {result.charge_range} includes charges:<br>{z_text}<br><br>"
            f"t_R: {result.tR:.2f} s<br>"
            f"σ: {result.sigma:.3e} s<br>"
            f"R²: {r2_colored} <br>"
            f"<b>R_h: {result.Rh:.3e} m</b><br>"
            f"D: {result.D:.3e} m²/s<br>"
            f"<a href='{tau_definition}' style='{plain_text_style}'>Tau</a>: {tau_colored}<br>"
            f"<a href='{peclet_definition}' style='{plain_text_style}'>Péclet</a>: {peclet_colored}<br>"
        )
        self.prev_btn.setEnabled(self.current_result_index > 0)
        self.next_btn.setEnabled(self.current_result_index < len(self.analysis_results) - 1)

    def is_row_valid(self, row: ProteinInputRow) -> bool:
        try:
            return all([
                row.mz.text().strip(),
                row.range.text().strip(),
                row.charge.text().strip(),
                row.charge_range.text().strip()
            ])
        except RuntimeError:
            return False

    def show_next(self):
        if self.current_result_index < len(self.analysis_results) - 1:
            self.current_result_index += 1
            self.show_current_result()

    def show_previous(self):
        if self.current_result_index > 0:
            self.current_result_index -= 1
            self.show_current_result()

    def calculate_tau_peclet(self):
        try:
            T = float(self.temp_input.text())
            viscosity = float(self.viscosity_input.text())
            R_h = float(self.estimated_rh_input.text())
            r_cap = float(self.radius_input.text())
            L = float(self.length_input.text())
            Q = float(self.flow_input.text())
        except ValueError:
            self.tau_out.setText("Invalid input ⚠️")
            self.peclet_out.setText("Invalid input ⚠️")
            return

        try:
            t = tau(T, L, viscosity, Q, R_h*10**-9)
            p = peclet(R_h*10**-9, T, viscosity, r_cap, Q)
        except Exception:
            self.tau_out.setText("Error")
            self.peclet_out.setText("Error")
            return

        self.tau_out.setText(f"{t:.4g}")
        self.peclet_out.setText(f"{p:.4g}")
        GREEN_BOX = "QLineEdit {"f"background-color: {accessibility_colors.GreenBox};border: 2px solid {accessibility_colors.GreenBox_Border};""color: black;font-weight: bold;}"

        Yellow_BOX = "QLineEdit {"f"background-color: {accessibility_colors.YellowBox};border: 2px solid {accessibility_colors.YellowBox_Border};color: black;""font-weight: bold;}"

        RED_BOX = "QLineEdit {"f"background-color: {accessibility_colors.RedBox};border: 2px solid {accessibility_colors.RedBox_Border};""color: black;font-weight: bold;}"

        from scipy.constants import k as boltzmann_c
        import math

        T_in_kelvin = T + 273.15
        L_in_meter = L * (10**-2)
        R_h_in_meter = R_h * (10**-9)
        radius_in_meter = r_cap * (10 ** -6)
        tau_threshold = 1.25
        Q_max_for_tau_m3s = (boltzmann_c * T_in_kelvin * L_in_meter) / (6 * viscosity * R_h_in_meter * tau_threshold)
        Q_max_for_tau_ulmin = (Q_max_for_tau_m3s * 1e9) * 60
        suggested_flow_rate_for_tau = math.floor(Q_max_for_tau_ulmin)

        peclet_threshold = 70
        Q_min_for_peclet_m3s = (peclet_threshold * boltzmann_c * T_in_kelvin * radius_in_meter) / (6 * viscosity * R_h_in_meter)
        Q_min_for_peclet_ulmin = (Q_min_for_peclet_m3s * 1e9) * 60
        suggested_flow_rate_for_peclet = math.ceil(Q_min_for_peclet_ulmin)

        if t > 1.25:
            self.tau_out.setText(self.tau_out.text()+"✔️")
            self.tau_out.setStyleSheet(GREEN_BOX)
        elif t>0.37:
            self.tau_out.setText(self.tau_out.text()+"⚠️")
            self.tau_out.setStyleSheet(Yellow_BOX)
            self.tau_out.setToolTip(
                f"Tau is below 1.25. Consider reducing the flow rate.\n"
                f"Suggested flow rate: ≤ {suggested_flow_rate_for_tau} µL/min to achieve τ ≥ 1.25."
            )
        else:
            self.tau_out.setText(self.tau_out.text()+"❌")
            self.tau_out.setStyleSheet(RED_BOX)
            self.tau_out.setToolTip(
                f"Tau is below 0.37. Reduce the flow rate significantly.\n"
                f"Suggested flow rate: ≤ {suggested_flow_rate_for_tau} µL/min to achieve τ ≥ 1.25."
            )

        if p > 70:
            self.peclet_out.setStyleSheet(GREEN_BOX)
            self.peclet_out.setText(self.peclet_out.text()+"✔️")
        else:
            self.peclet_out.setStyleSheet(RED_BOX)
            self.peclet_out.setText(self.peclet_out.text()+"❌")
            self.peclet_out.setToolTip(
                f"Péclet is below 70. Consider increasing the flow rate.\n"
                f"Suggested flow rate: ≥ {suggested_flow_rate_for_peclet} µL/min to achieve Péclet > 70."
            )

    def update_info(self, text: str):
        """
        Update the info panel on the right side of the UI.
        """
        self.info_box.setHtml(text)

    def export_pdf(self):
        if not self.analysis_results:
            QMessageBox.warning(
                self,
                "No results",
                "Run the analysis before exporting."
            )
            return

        now = datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        file_name = f"analysis_results_{timestamp}.pdf"

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export results as PDF",
            file_name,
            "PDF files (*.pdf)"
        )
        if not path:
            return

        try:
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.pyplot as plt
            import numpy as np

            temp = self.temp_input.text()
            viscosity = self.viscosity_input.text()
            radius = self.radius_input.text()
            length = self.length_input.text()
            flow = self.flow_input.text()
            #ms1_filename = os.path.basename(self.ms1_path)

            param_text = (
                f"File: {self.ms1_path}\n"
                f"T: {temp} °C  |  η: {viscosity} kg·m⁻¹·s⁻¹  |  "
                f"r: {radius} µm  |  L: {length} cm  |  Q: {flow} µL/min"
            )

            with PdfPages(path) as pdf:
                for i, result in enumerate(self.analysis_results, start=1):
                    fig = plt.figure(figsize=(8.5, 11))
                    gs = fig.add_gridspec(3, 1, height_ratios=[0.5,3, 1])

                    # ================= PARAMETERS HEADER =================
                    ax_header = fig.add_subplot(gs[0])
                    ax_header.axis("off")
                    ax_header.text(
                        0.5, 0.5, param_text,
                        ha="center", va="center",
                        fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="#e0e0e0", edgecolor="gray")
                    )

                    # ================= PLOT =================
                    ax_plot = fig.add_subplot(gs[1])

                    ax_plot.set_xlim(1, result.seconds[-1])
                    ax_plot.set_ylim(
                        min(*result.final_intensities, *result.removed_dip_fitted),
                        max(*result.final_intensities, *result.removed_dip_fitted)
                    )

                    # EIC
                    ax_plot.scatter(
                        result.seconds,
                        result.final_intensities,
                        s=20,
                        alpha=0.8,
                        color="red",
                        label="Masked data points"
                    )

                    # Masked
                    ax_plot.scatter(
                        result.seconds,
                        result.removed_dip,
                        s=20,
                        alpha=0.9,
                        color="orange",
                        label="EIC after Masking"
                    )

                    # Fit
                    ax_plot.plot(
                        result.seconds,
                        result.removed_dip_fitted,
                        "--",
                        linewidth=2,
                        color="blue",
                        label=f"Fit (R²={result.r2:.3f}) \nR_h = {(result.Rh * 10**9):.3f} nm "
                    )

                    ax_plot.set_xlabel("Time (s)")
                    ax_plot.set_ylabel("Intensity")
                    ax_plot.set_title(f"Protein {i}: {self.protein_rows[i-1].proteinName.text()} m/z {result.protein_mz} ",fontweight="bold")
                    ax_plot.grid(True)
                    ax_plot.legend(fontsize=9)


                    fig.subplots_adjust(
                        left=0.12,
                        right=0.95,
                        top=0.93,
                        bottom=0.32
                    )

                    # ================= TEXT =================
                    ax_text = fig.add_subplot(gs[2])
                    ax_text.axis("off")

                    summary = (
                        f"Protein {i} / {len(self.analysis_results)}\n"
                        f"{self.protein_rows[i-1].proteinName.text()}\n"
                        f"m/z: {result.protein_mz} range {result.mz_window}\n"
                        f"Charge: {result.charge_state} range {result.charge_range}\n\n"
                        f"t_R: {result.tR:.2f} s\n"
                        f"sigma: {result.sigma:.3e} s\n"
                        f"R²: {result.r2:.3f}\n"
                        f"R_h: {result.Rh:.3e} m\n"
                        f"D: {result.D} m²/s\n"
                        f"Tau: {result.t:.3f}\n"
                        f"Péclet: {result.p:.3e}"
                    )

                    ax_text.text(0.02, 0.95, summary, va="top", fontsize=11)

                    pdf.savefig(fig)
                    plt.close(fig)

        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def reset_for_accessibility(self):
        if hasattr(self,"analysis_results") and self.analysis_results != []:
            self.show_current_result()
            self.abort_remasking_btn.setEnabled(False)
            self.continue_remasking_btn.setEnabled(False)
        if hasattr(self,"tau_out") and self.tau_out.text() !="":
            self.calculate_tau_peclet()


class Accessibility_Window(QWidget):
    def __init__(self,parent,mode):
        super(Accessibility_Window,self).__init__(windowTitle="Accessibility",parent=parent)
        self.setWindowFlag(Qt.WindowType.Window)
        self.setWindowFlag(Qt.WindowType.MSWindowsFixedSizeDialogHint)
        layout = QVBoxLayout()
        self.title=QLabel("Accessibility Settings")
        self.drop_down=QComboBox()
        self.drop_down.addItems(["Normal","Deuteranomaly","Deuteranopia","Protanomaly","Protanopia","Tritanomaly","Tritanopia","Cone_Monochromacy","Achromatopsia"])
        self.drop_down.setCurrentText(mode)
        self.drop_down.currentIndexChanged.connect(self.on_change)
        layout.addWidget(self.title)
        layout.addWidget(self.drop_down)
        self.setLayout(layout)

    def on_change(self):
        accessibility_colors.change_accessibility_color(self.drop_down.currentText())
        self.parent().reset_for_accessibility()
