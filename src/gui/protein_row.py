from PySide6.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QPushButton

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
