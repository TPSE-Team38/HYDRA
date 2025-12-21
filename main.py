import sys
from PySide6.QtWidgets import QApplication
from src.gui.main_window import MainWindow

app = QApplication(sys.argv)
win = MainWindow()
win.show()
sys.exit(app.exec())
