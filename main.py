import multiprocessing
import sys

from PySide6.QtGui import QIcon,QBitmap
from PySide6.QtWidgets import QApplication
from src.gui.main_window import MainWindow

import sys

if __name__=="__main__":
    app = QApplication(sys.argv)
    icon=QIcon(QBitmap("./assets/HYDRA-logo.png"))
    app.setProperty("Main_Icon",icon)
    app.setWindowIcon(icon)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
