from matplotlib import figure
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar
)
from src.plotting import ResultPlot
from src.EIC_extraction import recalculate, results
import numpy as np
from src.models import AnalysisConfig
from src.pipeline import run_analysis


from src.Calculations import mz_to_mz


class PlotWidget(FigureCanvasQTAgg):
    def __init__(self):
        self.fig = Figure(figsize=(6, 4))
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)

    def show_eic(self, result, config= None):
        """
        result: object EICResult(seconds, final_intensities, removed_dip, removed_dip_fitted, r2, tR, sigma, D, Rh, t, p)
        """
        self.ax.clear()

        self.ax.set_xlim(1, result.seconds[-1])
        self.ax.set_ylim(min(result.final_intensities), max(result.final_intensities))
        self.ax.plot(result.seconds, result.final_intensities, lw=1.5)

        # Masking
        self.ax.scatter(result.seconds, result.removed_dip, label="EIC after Masking")
        self.ax.scatter(result.seconds, [v if v not in result.removed_dip else np.nan for v in result.final_intensities])

        #Fitting
        self.ax.set_ylim(min(*result.removed_dip_fitted, *result.final_intensities), max(*result.removed_dip_fitted, *result.final_intensities))
        self.ax.plot(result.seconds, result.removed_dip_fitted, '--',
                label=f"Fitted EIC with R^2 value of {result.r2} \n and R_h of {result.Rh}  \n D: {result.D} \n sigma: {result.sigma} \n tau: {result.t} \n peclet: {result.p}")

        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Intensity")
        self.ax.set_title("Extracted Ion Chromatogram")
        self.ax.grid(True)

        self.ax.legend(loc='best', fontsize=8)

        self.fig.tight_layout()
        self.draw()

        if config:
            params = [config.temperature, config.viscosity, config.capillary_radius, config.capillary_length, config.flow_rate, config.mz_window, config.charge_state, config.charge_range]

            results.append(ResultPlot(result.final_intensities, result.seconds, params, self.fig, self.ax, recalculate))

        #self.ax.canvas.draw()
        #self.ax.canvas.mpl_connect()