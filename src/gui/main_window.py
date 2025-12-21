from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QPushButton,
    QFileDialog, QLabel, QTextEdit,
    QLineEdit, QFormLayout, QSplitter, QMessageBox
)
from PySide6.QtCore import Qt

from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from .plot_widget import PlotWidget
from .controller import AnalysisController
from src.models import AnalysisConfig


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Catalyst – EIC Analysis")
        self.resize(1100, 700)

        self.ms1_path = None

        # Use PlotWidget directly (it already contains Figure and Canvas!)
        self.plot = PlotWidget()

        # Initialize controller with the PlotWidget
        self.controller = AnalysisController(self.plot, self)

        self._build_ui()

    def _build_ui(self):
        main_layout = QVBoxLayout()

        # File loader
        file_layout = QHBoxLayout()
        self.file_label = QLabel("No MS1 file selected")

        load_btn = QPushButton("Load MS1 file")
        load_btn.clicked.connect(self.load_file)

        file_layout.addWidget(load_btn)
        file_layout.addWidget(self.file_label)

        # Parameters
        form = QFormLayout()

        self.mz_input = QLineEdit("1689")
        self.window_input = QLineEdit("2")
        self.charge_input = QLineEdit("7")
        self.charge_range_input = QLineEdit("2")

        self.temp_input = QLineEdit("298")
        self.viscosity_input = QLineEdit("0.001")
        self.radius_input = QLineEdit("75e-6")
        self.length_input = QLineEdit("0.5")
        self.flow_input = QLineEdit("20e-9")

        form.addRow("Protein m/z", self.mz_input)
        form.addRow("m/z window (±)", self.window_input)
        form.addRow("Charge state", self.charge_input)
        form.addRow("Charge range", self.charge_range_input)

        form.addRow("Temperature (K)", self.temp_input)
        form.addRow("Viscosity (Pa·s)", self.viscosity_input)
        form.addRow("Capillary radius (m)", self.radius_input)
        form.addRow("Capillary length (m)", self.length_input)
        form.addRow("Flow rate (m³/s)", self.flow_input)

        # Run button
        run_btn = QPushButton("Run analysis")
        run_btn.clicked.connect(self.run)

        # Splitter: plot | info
        splitter = QSplitter(Qt.Horizontal)

        # Plot side - create container widget
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)

        # Add toolbar for zoom/pan functionality
        # PlotWidget is already a FigureCanvas, so just pass it to the toolbar
        self.toolbar = NavigationToolbar(self.plot, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.plot)  # self.plot is the PlotWidget (canvas)

        # Info side
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setMaximumWidth(350)

        # Add both sides to splitter
        splitter.addWidget(plot_widget)
        splitter.addWidget(self.info_box)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Assemble main layout - add splitter ONCE at the end
        main_layout.addLayout(file_layout)
        main_layout.addLayout(form)
        main_layout.addWidget(run_btn)
        main_layout.addWidget(splitter)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        self.update_info("Load an MS1 file to begin.")

    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select MS1 file",
            "",
            "MS1 files (*.ms1)"
        )

        if path:
            self.ms1_path = path
            self.file_label.setText(path.split("/")[-1])

    def run(self):
        if not self.ms1_path:
            QMessageBox.warning(self, "Error", "Please load an MS1 file first.")
            return

        try:
            config = AnalysisConfig(
                ms1_path=self.ms1_path,
                protein_mz=float(self.mz_input.text()),
                mz_window=float(self.window_input.text()),
                charge_state=int(self.charge_input.text()),
                charge_range=int(self.charge_range_input.text()),
                temperature=float(self.temp_input.text()),
                viscosity=float(self.viscosity_input.text()),
                capillary_radius=float(self.radius_input.text()),
                capillary_length=float(self.length_input.text()),
                flow_rate=float(self.flow_input.text())
            )

            self.controller.run(config)

        except ValueError:
            QMessageBox.critical(self, "Invalid input", "Please check parameter values.")

    def update_info(self, text):
        self.info_box.setText(text)