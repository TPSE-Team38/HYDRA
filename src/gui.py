import sys
import numpy as np
from pyteomics import ms1
#from EIC_extraction import get_final_eic_intensities
#from EIC_extraction import load_ms1

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QPushButton,
    QFileDialog, QLabel, QTextEdit,
    QLineEdit, QFormLayout, QSplitter
)
from PySide6.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar
)



# Backend: Minimal EIC extraction
# we will use the eic_extraction.py code in the future

class EICExtractor:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.spectra = []

    def load_ms1(self):
        """Load MS1 spectra"""
        self.spectra.clear()
        with open(self.filepath, "r") as f:
            for spectrum in ms1.read(f):
                self.spectra.append(spectrum)

    def extract_eic(self, target_mz: float, mz_window: float):
        """
        Extract EIC for ONE m/z region
        Returns:
            x : seconds
            y : summed intensity
        """
        intensities = []

        for spectrum in self.spectra:
            mz = spectrum["m/z array"]
            inten = spectrum["intensity array"]

            mask = (mz >= target_mz - mz_window) & (mz <= target_mz + mz_window)
            intensities.append(np.sum(inten[mask]))

        x = np.arange(1, len(intensities) + 1)
        y = np.array(intensities)

        return x, y


# =========================================================
# GUI
# =========================================================
class EICExtractorGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Minimal EIC Extractor (MS1)")
        self.setGeometry(100, 100, 1200, 800)

        self.current_file = None
        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # -------------------------------
        # Top buttons
        # -------------------------------
        top_layout = QHBoxLayout()

        self.load_btn = QPushButton("Load MS1 file")
        self.load_btn.clicked.connect(self.load_file)

        self.run_btn = QPushButton("Extract EIC")
        self.run_btn.clicked.connect(self.run_eic)
        self.run_btn.setEnabled(False)

        self.file_label = QLabel("No file loaded")

        top_layout.addWidget(self.load_btn)
        top_layout.addWidget(self.run_btn)
        top_layout.addWidget(self.file_label)
        top_layout.addStretch()

        main_layout.addLayout(top_layout)

        # -------------------------------
        # Parameters
        # -------------------------------
        param_layout = QFormLayout()

        self.mz_input = QLineEdit("1689")
        self.range_input = QLineEdit("2")

        param_layout.addRow("Target m/z:", self.mz_input)
        param_layout.addRow("m/z window (±):", self.range_input)

        main_layout.addLayout(param_layout)

        # -------------------------------
        # Splitter: plot | info
        # -------------------------------
        splitter = QSplitter(Qt.Horizontal)

        # Plot side
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)

        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        # Info side
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setMaximumWidth(350)

        splitter.addWidget(plot_widget)
        splitter.addWidget(self.info_box)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        main_layout.addWidget(splitter)

        self.update_info("Load an MS1 file to begin.")

    # -------------------------------------------------
    # Actions
    # -------------------------------------------------
    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select MS1 file",
            "",
            "MS1 files (*.ms1)"
        )

        if file_path:
            self.current_file = file_path
            self.file_label.setText(file_path.split("/")[-1])
            self.run_btn.setEnabled(True)
            self.update_info(f"Loaded file:\n{file_path}")

    def run_eic(self):
        try:
            target_mz = float(self.mz_input.text())
            mz_window = float(self.range_input.text())

            extractor = EICExtractor(self.current_file)
            extractor.load_ms1()
            self.spe
            extractor.load_ms1()
            x, y = extractor.extract_eic(target_mz, mz_window)


            self.plot_eic(x, y)

            self.update_info(
                f"EIC extracted successfully\n\n"
                f"Target m/z: {target_mz}\n"
                f"Window: ±{mz_window}\n"
                f"Scans: {len(x)}\n"
                f"Max intensity: {y.max():.2e}"
            )

        except Exception as e:
            self.update_info(f"Error:\n{e}")

    # -------------------------------------------------
    # Plotting
    # -------------------------------------------------
    def plot_eic(self, x, y):
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        ax.plot(x, y, lw=1.5)
        ax.set_xlabel("Seconds")
        ax.set_ylabel("Summed intensity")
        ax.set_title("Extracted Ion Chromatogram (EIC)")
        ax.grid(True)


        self.figure.tight_layout()
        self.canvas.draw()
        self.canvas.mpl_connect()

    def update_info(self, text):
        self.info_box.setText(text)


# =========================================================
# Entry point
# =========================================================
def main():
    app = QApplication(sys.argv)
    window = EICExtractorGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
