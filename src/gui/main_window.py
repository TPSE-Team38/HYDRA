from datetime import datetime
from pathlib import Path

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QLabel, QPushButton, QFileDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox,
    QLineEdit, QSplitter, QTextEdit, QMessageBox, QScrollArea, QApplication,QSystemTrayIcon
)
from PySide6.QtCore import Qt,QThreadPool,Signal,QObject,QRunnable
from PySide6.QtGui import QIcon,QPixmap

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from .plot_widget import PlotWidget
from .controller import AnalysisController
from .protein_row import ProteinInputRow
from src.models import AnalysisConfig
from src.Calculations import tau, peclet

import os

from ..parallization import LoadMS1Worker


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Catalyst 2.0")
        self.resize(1200, 800)
        self.reset_btn,self.abort_remasking_btn,self.continue_remasking_btn =None,None,None
        self.icon=QIcon("./assets/catalyst.ico")
        self.setWindowIcon(self.icon)
        self.ms1_path = None
        self.protein_rows = []
        self.analysis_results = []  # list of results
        self.current_result_index = 0  # which protein is shown

        self.plot = PlotWidget()
        self.controller = AnalysisController(self.plot, self)
        self._build_ui()

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

        title = QLabel("CATALYST 2.0")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 21px; font-weight: bold;")

        top_bar.addWidget(self.load_btn)
        top_bar.addWidget(self.file_status)
        top_bar.addStretch()
        top_bar.addWidget(title)
        top_bar.addStretch()

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
        plot_layout = QVBoxLayout(plot_container)

        # Remove margins so the plot fits snugly
        plot_layout.setContentsMargins(0, 0, 0, 0)

        toolbar = NavigationToolbar(self.plot, self)

        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Results Plot"))

        open_fullscreen_btn = QPushButton("View Plot in Full Screen")
        open_fullscreen_btn.clicked.connect(self.open_plot_fullscreen)

        header_layout.addStretch()
        header_layout.addWidget(open_fullscreen_btn)

        plot_layout.addLayout(header_layout)
        plot_layout.addWidget(toolbar)
        plot_layout.addWidget(self.plot)

        lower_splitter.addWidget(plot_container)

        # main_splitter.setStretchFactor(0, 0)
        # main_splitter.setStretchFactor(1, 1)

        # info box
        self.info_box = QTextEdit()
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
        self.continue_remasking_btn.setEnabled(False)
        self.abort_remasking_btn.setEnabled(False)
        # ---------- Reset Button ----------
        self.reset_btn = QPushButton("Reset Masking")
        self.reset_btn.setEnabled(True)
        nav_layout.addWidget(self.prev_btn)
        # nav_layout.addStretch()
        nav_layout.addWidget(self.abort_remasking_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.reset_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.continue_remasking_btn)
        # nav_layout.addStretch(1)
        nav_layout.addWidget(self.next_btn)

        main_layout.addLayout(nav_layout)


    # ================= ACTIONS =================
    def open_plot_fullscreen(self):
        if self.plot.stored_show_eic_args is None:
            QMessageBox.warning(self, "No plot", "No plot available yet.")
            return

        self.plot_window = QMainWindow(self)
        self.plot_window.setWindowTitle("Results Plot")
        self.plot_window.setWindowIcon(self.icon)
        self.plot_window.resize(1400, 900)

        central = QWidget()
        layout = QVBoxLayout(central)

        # new PlotWidget (same class)
        fullscreen_plot = PlotWidget()

        toolbar = NavigationToolbar(fullscreen_plot, self.plot_window)

        layout.addWidget(toolbar)
        layout.addWidget(fullscreen_plot)

        self.plot_window.setCentralWidget(central)

        # redraw using same data
        fullscreen_plot.show_eic(*self.plot.stored_show_eic_args)

        self.plot_window.showMaximized()

    def reset_masking(self):
        self.show_current_result()

    def show_recalculated_fit(self,masked_y, fitted_y, r2, t_R, sigma, D, R_h, t, p):
        result=self.analysis_results[self.current_result_index]
        self.update_info(
            f"Protein {self.current_result_index + 1} / {len(self.analysis_results)}\n"
            f"{self.protein_rows[self.current_result_index].proteinName.text()}\n"
            f"m/z: {result.protein_mz} range {result.mz_window}\n"
            f"Charge: {result.charge_state} range {result.charge_range}\n\n"
            "Previous Fit:\n \n"
            f"t_R: {result.tR:.2f} s\n"
            f"σ: {result.sigma:.3e}\n"
            f"R²: {result.r2: }\n"
            f"R_h: {result.Rh:.3e} m\n"
            f"D: {result.D:.3e} m²/s\n"
            f"Tau: {result.t:.3f}\n"
            f"Péclet: {result.p:.3e}\n \n"
            "Recalculated Fit:\n"
            f"t_R: {t_R:.2f} s\n"
            f"σ: {sigma:.3e}\n"
            f"R²: {r2:}\n"
            f"R_h: {R_h:.3e} m\n"
            f"D: {D:.3e} m²/s\n"
            f"Tau: {t:.3f}\n"
            f"Péclet: {p:.3e}\n"
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
        self.file_status.setText("✔ ms1 file loaded")
        self.file_status.setStyleSheet("color: green;")

    def on_ms1_error(self, msg):
        self.file_status.setText("❌ Failed to load ms1")
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
        self.plot.show_eic(result,self.reset_btn,self.show_current_result,self.abort_remasking_btn,self.continue_remasking_btn,self.show_recalculated_fit,config=AnalysisConfig(
            ms1_path=self.ms1_path,
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

        self.update_info(
            f"Protein {self.current_result_index + 1} / {len(self.analysis_results)}\n"
            f"{self.protein_rows[self.current_result_index].proteinName.text()}\n"
            f"m/z: {result.protein_mz} range {result.mz_window}\n"
            f"Charge: {result.charge_state} range {result.charge_range}\n\n"
            f"t_R: {result.tR:.2f} s\n"
            f"σ: {result.sigma:.3e}\n"
            f"R²: {result.r2: }\n"
            f"R_h: {result.Rh:.3e} m\n"
            f"D: {result.D:.3e} m²/s\n"
            f"Tau: {result.t:.3f}\n"
            f"Péclet: {result.p:.3e}"
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
            self.tau_out.setText("Invalid input")
            self.peclet_out.setText("Invalid input")
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
        GREEN_BOX = """
        QLineEdit {
            background-color: #c6f6d5;
            border: 2px solid #2e7d32;
            color: black;
            font-weight: bold;
        }
        """
        RED_BOX = """
        QLineEdit {
            background-color: #fed7d7;
            border: 2px solid #c62828;
            color: black;
            font-weight: bold;
        }
        """
        if t > 1.25:
            self.tau_out.setStyleSheet(GREEN_BOX)
        else:
            self.tau_out.setStyleSheet(RED_BOX)

        if p > 70:
            self.peclet_out.setStyleSheet(GREEN_BOX)
        else:
            self.peclet_out.setStyleSheet(RED_BOX)

    def update_info(self, text: str):
        """
        Update the info panel on the right side of the UI.
        """
        self.info_box.setText(text)

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

            with PdfPages(path) as pdf:
                for i, result in enumerate(self.analysis_results, start=1):
                    fig = plt.figure(figsize=(8.5, 11))
                    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])

                    # ================= PLOT =================
                    ax_plot = fig.add_subplot(gs[0])

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
                        label=f"Fit (R²={result.r2:.3f})"
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
                    ax_text = fig.add_subplot(gs[1])
                    ax_text.axis("off")

                    summary = (
                        f"Protein {i} / {len(self.analysis_results)}\n"
                        f"{self.protein_rows[i-1].proteinName.text()}\n"
                        f"m/z: {result.protein_mz} range {result.mz_window}\n"
                        f"Charge: {result.charge_state} range {result.charge_range}\n\n"
                        f"t_R: {result.tR:.2f} s\n"
                        f"sigma: {result.sigma:.3e}\n"
                        f"R²: {result.r2:.3f}\n"
                        f"R_h: {result.Rh:.3e} m\n"
                        f"Tau: {result.t:.3f}\n"
                        f"Péclet: {result.p:.3e}"
                    )

                    ax_text.text(0.02, 0.95, summary, va="top", fontsize=11)

                    pdf.savefig(fig)
                    plt.close(fig)

        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

