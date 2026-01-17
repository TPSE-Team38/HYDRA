import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication
from src.gui.main_window import MainWindow

app = QApplication(sys.argv)
icon=QIcon()
icon.addFile("./assets/catalyst.ico")
app.setWindowIcon(icon)
win = MainWindow()
win.show()
sys.exit(app.exec())
