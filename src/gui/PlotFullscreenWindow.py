from PySide6.QtWidgets import QMainWindow, QMessageBox


class PlotFullscreenWindow(QMainWindow):
    def __init__(self, parent, plot, toolbar, on_close):
        super().__init__(parent)
        self.plot = plot
        self.toolbar = toolbar
        self.on_close = on_close
        self.parent_window = parent

    def closeEvent(self, event):
        # If there is an active remasking session, block closing
        if (
                self.parent_window.abort_remasking_btn.isEnabled()
                or self.parent_window.continue_remasking_btn.isEnabled()):
            QMessageBox.warning(
                self,
                "Confirm or Abort Remasking",
                "You must confirm or abort remasking before closing the window."
            )
            event.ignore()  # block closing

        else:
            self.on_close()
            event.accept()


