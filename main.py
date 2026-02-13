import multiprocessing
import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication
from src.gui.main_window import MainWindow

import sys

# Add this BEFORE freeze_support()
# if sys.platform.startswith('win'):
#     # Workaround for PyInstaller + multiprocessing
#     multiprocessing.set_start_method('spawn', force=True)
if __name__=="__main__":
    # multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    icon=QIcon()
    icon.addFile("./assets/catalyst.ico")
    app.setProperty("Main_Icon",icon)
    app.setWindowIcon(icon)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
