from PySide6.QtWidgets import QMainWindow

class PlotFullscreenWindow(QMainWindow):
    def __init__(self, parent, plot, toolbar, on_close):
        super().__init__(parent)
        self.plot = plot
        self.toolbar = toolbar
        self.on_close = on_close

    def closeEvent(self, event):
        self.on_close()
        event.accept()
