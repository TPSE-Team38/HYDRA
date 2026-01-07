from src.pipeline import run_analysis

class AnalysisController:
    def __init__(self, plot, window):
        self.plot = plot
        self.window = window

    def run(self, config, store=True):
        result = run_analysis(config)  # backend

        if store:
            self.plot.show_eic(result)

        return result
