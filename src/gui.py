from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QLabel, QPushButton, QFileDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLineEdit, QSplitter, QScrollArea
)
from PySide6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from src.EIC_extraction import load_ms1
from src.ms1_analysis import analyse
from src.Calculations import tau, peclet


# ================= PROTEIN INPUT ROW =================
class ProteinInputRow(QWidget):
    def __init__(self):
        super().__init__()

        layout = QHBoxLayout(self)

        self.mz = QLineEdit()
        self.mz.setPlaceholderText("Protein m/z")

        self.range = QLineEdit()
        self.range.setPlaceholderText("Sampling range")

        self.charge = QLineEdit()
        self.charge.setPlaceholderText("Charge state")

        self.charge_range = QLineEdit()
        self.charge_range.setPlaceholderText("Charge range")

        remove_btn = QPushButton("✕")
        remove_btn.setFixedWidth(28)
        remove_btn.clicked.connect(self.remove_self)

        layout.addWidget(self.mz)
        layout.addWidget(self.range)
        layout.addWidget(self.charge)
        layout.addWidget(self.charge_range)
        layout.addWidget(remove_btn)

        self.analysis_results = []
        self.current_result_index = 0

    def remove_self(self):
        self.setParent(None)
        self.deleteLater()


# ================= GRAPH WINDOW =================
class GraphWindow(QMainWindow):
    def __init__(self, graph_widget):
        super().__init__()
        self.setWindowTitle("Results Plot")
        self.setCentralWidget(graph_widget)
        self.resize(1000, 700)


# ================= MAIN WINDOW =================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Catalyst 2.0")
        self.setMinimumSize(900, 600)
        self.showMaximized()

        self.spectra = None
        self.graph_window = None
        self.protein_rows = []
        self.analysis_results = None
        self.canvas = None

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(8, 8, 8, 8)


        # ================= TOP BAR =================
        top_bar = QWidget()
        top_layout = QGridLayout(top_bar)
        top_layout.setContentsMargins(10, 5, 10, 5)

        left_widget = QWidget()
        left_layout = QHBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)

        upload_btn = QPushButton("Upload ms1 file")
        upload_btn.setMinimumSize(140, 40)
        upload_btn.clicked.connect(self.upload_ms1)

        self.ms1_status = QLabel("No file loaded")
        self.ms1_status.setStyleSheet("color: red;")

        left_layout.addWidget(upload_btn)
        left_layout.addWidget(self.ms1_status)

        title = QLabel("CATALYST 2.0")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 21px; font-weight: bold;")

        right_spacer = QWidget()
        right_spacer.setFixedWidth(left_widget.sizeHint().width())

        top_layout.addWidget(left_widget, 0, 0)
        top_layout.addWidget(title, 0, 1)
        top_layout.addWidget(right_spacer, 0, 2)

        top_layout.setColumnStretch(0, 0)
        top_layout.setColumnStretch(1, 1)
        top_layout.setColumnStretch(2, 0)

        main_layout.addWidget(top_bar)

        # ================= MAIN SPLITTER =================
        self.main_splitter = QSplitter(Qt.Orientation.Vertical)
        main_layout.addWidget(self.main_splitter)

        # ================= UPPER HALF =================
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

        self.temp_input = QLineEdit()
        self.viscosity_input = QLineEdit()
        self.capillary_radius_input = QLineEdit()
        self.capillary_length_input = QLineEdit()
        self.flow_rate_input = QLineEdit()
        self.estimated_rh_input = QLineEdit()

        labels_inputs = [
            ("Temperature", self.temp_input),
            ("Viscosity", self.viscosity_input),
            ("Capillary radius", self.capillary_radius_input),
            ("Capillary length", self.capillary_length_input),
            ("Flow rate", self.flow_rate_input),
            ("Estimated Rh", self.estimated_rh_input)
        ]

        for row, (label, widget) in enumerate(labels_inputs):
            param_layout.addWidget(QLabel(label + ":"), row, 0)
            param_layout.addWidget(widget, row, 1)


        # ---------- Tau & Péclet ----------
        tau_box = QGroupBox()
        tau_layout = QVBoxLayout(tau_box)

        self.calc_tau_btn = QPushButton("Calculate Tau and Péclet")
        self.calc_tau_btn.setStyleSheet("font-weight: bold;")
        self.calc_tau_btn.clicked.connect(self.calculate_tau_peclet)

        tau_layout.addWidget(self.calc_tau_btn)

        tau_layout.addWidget(QLabel("Tau(τ):"))
        self.tau_output = QLineEdit()
        self.tau_output.setReadOnly(True)
        tau_layout.addWidget(self.tau_output)

        tau_layout.addWidget(QLabel("Péclet:"))
        self.peclet_output = QLineEdit()
        self.peclet_output.setReadOnly(True)
        tau_layout.addWidget(self.peclet_output)

        tau_layout.addStretch()

        upper_splitter.addWidget(protein_box)
        upper_splitter.addWidget(param_box)
        upper_splitter.addWidget(tau_box)

        upper_splitter.setStretchFactor(0, 1)
        upper_splitter.setStretchFactor(1, 1)
        upper_splitter.setStretchFactor(2, 1)

        upper_layout.addWidget(upper_splitter, stretch=1)


        # ---------- Analyse Button ----------

        analyse_btn = QPushButton("Analyse Data")
        analyse_btn.setMinimumSize(130, 30)
        analyse_btn.setStyleSheet("font-size: 16px; font-weight: bold;")
        analyse_btn.clicked.connect(self.on_analyse_clicked)

        analyse_row = QHBoxLayout()
        analyse_row.addStretch()
        analyse_row.addWidget(analyse_btn)
        analyse_row.addStretch()

        upper_layout.addLayout(analyse_row)

        # ================= LOWER HALF =================
        self.graph_container = QWidget()
        self.graph_container.setMinimumHeight(300)
        self.graph_layout = QVBoxLayout(self.graph_container)

        graph_header = QHBoxLayout()
        graph_title = QLabel("Results Plot")

        popout_btn = QPushButton("⬈")
        popout_btn.setFixedSize(32, 32)
        popout_btn.clicked.connect(self.popout_graph)

        graph_header.addWidget(graph_title)
        graph_header.addStretch()
        graph_header.addWidget(popout_btn)

        self.graph_layout.addLayout(graph_header)

        self.graph_placeholder = QLabel("Graph will appear here")
        self.graph_placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.graph_placeholder.setStyleSheet("border: 1px dashed gray;")

        self.graph_layout.addWidget(self.graph_placeholder, stretch=1)

        self.main_splitter.addWidget(upper_widget)
        self.main_splitter.addWidget(self.graph_container)
        self.main_splitter.setStretchFactor(0, 0)
        self.main_splitter.setStretchFactor(1, 1)

        # --- Navigation buttons ---
        nav_layout = QHBoxLayout()

        self.prev_btn = QPushButton("< Previous")
        self.next_btn = QPushButton("Next >")

        self.prev_btn.clicked.connect(self.show_previous_plot)
        self.next_btn.clicked.connect(self.show_next_plot)

        nav_layout.addWidget(self.prev_btn)
        nav_layout.addStretch()
        nav_layout.addWidget(self.next_btn)

        self.graph_layout.addLayout(nav_layout)


    # ================= FUNCTIONS =================

    def add_protein_row(self):
        row = ProteinInputRow()
        self.protein_rows.append(row)
        self.protein_list_layout.addWidget(row)

    def upload_ms1(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select ms1 file", "", "ms1 files (*.ms1)"
        )
        if not path:
            return

        self.ms1_status.setText("Loading...")
        self.ms1_status.setStyleSheet("color: orange; font-weight: bold;")
        QApplication.processEvents()

        self.spectra = load_ms1(path)

        self.ms1_status.setText("✔ ms1 file loaded")
        self.ms1_status.setStyleSheet("color: green; font-weight: bold;")

    def popout_graph(self):
        if self.graph_window is None:
            self.graph_window = GraphWindow(self.graph_container)
            self.graph_window.showMaximized()

    def on_analyse_clicked(self):
        if not hasattr(self, "spectra") or self.spectra is None:
            print("No ms1 file loaded")
            return

        try:
            parameters = {
                "temperature": float(self.temp_input.text()),
                "viscosity": float(self.viscosity_input.text()),
                "capillary_radius": float(self.capillary_radius_input.text()),
                "capillary_length": float(self.capillary_length_input.text()),
                "flow_rate": float(self.flow_rate_input.text())
            }
        except ValueError:
            print("Invalid parameter input")
            return

        # ---------- Read protein inputs ----------
        regions = []

        for i, row in enumerate(self.protein_rows, start=1):

            try:
                regions.append((
                    float(row.mz.text()),
                    float(row.range.text()),
                    int(row.charge.text()),
                    int(row.charge_range.text()),
                    f"Protein {i}"
                ))
            except ValueError:
                print(f"Invalid protein input in row {i}")
                return

        if not regions:
            print("No valid proteins entered")
            return

        # ---------- Run analysis ----------
        results = analyse(
            spectra=self.spectra,
            parameters=parameters,
            regions=regions
        )

        self.analysis_results = results
        self.current_result_index = 0

        # Show first protein
        self.show_figure(self.analysis_results[0].fig)
        self.main_splitter.setSizes([210, 790])#???????

    def show_next_plot(self):
        if not self.analysis_results:
            return

        if self.current_result_index < len(self.analysis_results) - 1:
            self.current_result_index += 1
            self.show_figure(self.analysis_results[self.current_result_index].fig)

    def show_previous_plot(self):
        if not self.analysis_results:
            return

        if self.current_result_index > 0:
            self.current_result_index -= 1
            self.show_figure(self.analysis_results[self.current_result_index].fig)

    def show_figure(self, fig):
        # Remove placeholder
        self.graph_placeholder.hide()

        # Clear old canvas
        if hasattr(self, "canvas") and self.canvas is not None:
            self.graph_layout.removeWidget(self.canvas)
            self.canvas.setParent(None)

        self.canvas = FigureCanvas(fig)
        self.graph_layout.insertWidget(1, self.canvas, stretch=1)

    def calculate_tau_peclet(self):
        try:
            T = float(self.temp_input.text())
            viscosity = float(self.viscosity_input.text())
            R_h = float(self.estimated_rh_input.text())
            r_cap = float(self.capillary_radius_input.text())
            L = float(self.capillary_length_input.text())
            Q = float(self.flow_rate_input.text())
        except ValueError:
            self.tau_output.setText("Invalid input")
            self.peclet_output.setText("Invalid input")
            return

        try:
            t = tau(T, L, viscosity, Q, R_h)
            p = peclet(R_h, T, viscosity, r_cap, Q)
        except Exception:
            self.tau_output.setText("Error")
            self.peclet_output.setText("Error")
            return

        self.tau_output.setText(f"{t:.4g}")
        self.peclet_output.setText(f"{p:.4g}")
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
            self.tau_output.setStyleSheet(GREEN_BOX)
        else:
            self.tau_output.setStyleSheet(RED_BOX)

        if p > 70:
            self.peclet_output.setStyleSheet(GREEN_BOX)
        else:
            self.peclet_output.setStyleSheet(RED_BOX)


