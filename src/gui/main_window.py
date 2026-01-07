from PySide6.QtWidgets import (
    QMainWindow, QWidget, QLabel, QPushButton, QFileDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox,
    QLineEdit, QSplitter, QTextEdit, QMessageBox, QScrollArea
)
from PySide6.QtCore import Qt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from .plot_widget import PlotWidget
from .controller import AnalysisController
from .protein_row import ProteinInputRow
from src.models import AnalysisConfig
from src.Calculations import tau, peclet


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Catalyst 2.0")
        self.resize(1200, 800)

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

        load_btn = QPushButton("Upload ms1 file")
        load_btn.clicked.connect(self.load_file)

        self.file_status = QLabel("No ms1 file loaded")
        self.file_status.setStyleSheet("color: red;")

        title = QLabel("CATALYST 2.0")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 21px; font-weight: bold;")

        top_bar.addWidget(load_btn)
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
        analyse_btn = QPushButton("Analyse Data")
        analyse_btn.setStyleSheet("font-size: 16px; font-weight: bold;")
        analyse_btn.clicked.connect(self.run)

        upper_layout.addWidget(analyse_btn, alignment=Qt.AlignmentFlag.AlignCenter)

        main_splitter.addWidget(upper_widget)

        # ================= LOWER PANEL =================
        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)

        toolbar = NavigationToolbar(self.plot, self)

        plot_layout.addWidget(QLabel("Results Plot"))
        plot_layout.addWidget(toolbar)
        plot_layout.addWidget(self.plot)

        main_splitter.addWidget(plot_container)

        main_splitter.setStretchFactor(0, 0)
        main_splitter.setStretchFactor(1, 1)

        self.setCentralWidget(central)

        # Navigation buttons

        nav_layout = QHBoxLayout()

        self.prev_btn = QPushButton("< Previous")
        self.next_btn = QPushButton("Next >")

        self.prev_btn.clicked.connect(self.show_previous)
        self.next_btn.clicked.connect(self.show_next)

        nav_layout.addWidget(self.prev_btn)
        nav_layout.addStretch()
        nav_layout.addWidget(self.next_btn)

        main_layout.addLayout(nav_layout)

        # info box
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setMinimumWidth(250)
        self.info_box.setPlaceholderText("Analysis results will appear here")
        main_splitter.addWidget(self.info_box)
        main_splitter.setStretchFactor(1, 1)

    # ================= ACTIONS =================

    def add_protein_row(self):
        row = ProteinInputRow()
        self.protein_rows.append(row)
        self.protein_list_layout.addWidget(row)

    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select ms1 file", "", "*.ms1")
        if path:
            self.ms1_path = path
            self.file_status.setText("✔ ms1 file loaded")
            self.file_status.setStyleSheet("color: green;")

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

                result = self.controller.run(config, store=False)
                if result is None:
                    QMessageBox.critical(self, "Analysis failed", "Analysis returned no result.")
                    return
                self.analysis_results.append(result)

            self.current_result_index = 0
            self.show_current_result()

        except ValueError:
            QMessageBox.critical(self, "Invalid input", "Check protein or parameter values.")

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

    def show_current_result(self):
        if not self.analysis_results:
            return

        result = self.analysis_results[self.current_result_index]
        self.plot.show_eic(result)

        self.update_info(
            f"Protein {self.current_result_index + 1} / {len(self.analysis_results)}\n"
            f"m/z: {result.protein_mz} range {result.mz_window}\n"
            f"Charge: {result.charge_state} range {result.charge_range}\n\n"
            f"t_R: {result.tR:.2f} s\n"
            f"σ: {result.sigma:.3e}\n"
            f"r2: {result.r2: }\n"
            f"Rh: {result.Rh:.3e} m\n"
            f"Tau: {result.t:.3f}\n"
            f"Péclet: {result.p:.3e}"
        )
        self.prev_btn.setEnabled(self.current_result_index > 0)
        self.next_btn.setEnabled(self.current_result_index < len(self.analysis_results) - 1)

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
            t = tau(T, L, viscosity, Q, R_h)
            p = peclet(R_h, T, viscosity, r_cap, Q)
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
