from typing import Callable

from PySide6.QtCore import Signal,QObject,QRunnable,Slot

class WorkerSignals(QObject):
    finished = Signal(str)
    error = Signal(str)

class LoadMS1Worker(QRunnable):
    def __init__(self, func:Callable, path):
        super().__init__()
        self.func = func
        self.path = path
        self.signals = WorkerSignals()

    @Slot()
    def run(self):
        try:
            self.func(self.path)
            self.signals.finished.emit(self.path)
        except Exception as e:
            self.signals.error.emit(str(e))