from src.pipeline import run_analysis

class AnalysisController:
    def __init__(self, plot_widget, parent):
        self.plot = plot_widget
        self.parent = parent

    def run(self, config):

        result = run_analysis(config)
        self.plot.show_eic(result)

