from src.EIC_extraction import load_ms1, load_ms1_parallel
from src.pipeline import run_analysis

class AnalysisController:
    def __init__(self, plot, window):
        self.plot = plot
        self.window = window
        self.spectra = None
        self.ms1_path = None

    def load_ms1_once(self, path):
        if self.spectra is None or self.ms1_path != path:
            self.spectra = load_ms1_parallel(path)
            self.ms1_path = path

    def run(self, config,reset_btn,abort_remasking_btn,continue_remasking_btn,show_result, store=True):
        if self.spectra is None:
            raise RuntimeError("ms1 file not loaded")

        result = run_analysis(self.spectra, config)

        # if store:
        self.plot.show_eic(result,reset_btn,abort_remasking_btn,continue_remasking_btn,show_result, config)

        return result
