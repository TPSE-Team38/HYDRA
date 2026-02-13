from typing import Callable

import matplotlib.backend_bases
import matplotlib.pyplot as plt
import numpy as np
from PySide6.QtWidgets import QPushButton

from matplotlib.lines import Line2D
from matplotlib.widgets import Button
from src.models import EICResult

class ResultPlot:
    y: np.ndarray
    x: np.ndarray
    peaks: list[tuple[float,float]]
    recalculated_fit: Line2D
    params: dict
    fig: plt.Figure
    ax: plt.Axes
    func: Callable
    new_points_plot: list[plt.Line2D]

    def __init__(self,y,x,params,fig,ax,func:Callable,EIC_result:EICResult,reset_btn:QPushButton,abort_remasking_btn:QPushButton,continue_remasking_btn:QPushButton,show_result:Callable,show_recalculated_fit:Callable):
        self.results_of_remasking = None
        self.continue_btn = continue_remasking_btn
        self.abort_btn = abort_remasking_btn
        self.peaks=[]
        self.y=y
        self.x=x
        self.new_points_plot=[]
        self.params=params
        self.fig,self.ax=fig,ax
        self.cid=self.fig.canvas.mpl_connect("button_press_event",self.on_click)
        self.recalculated_fit = None
        self.func=func
        self.EIC_result = EIC_result
        self.reset_btn=reset_btn
        self.reset_btn.clicked.connect(self.on_reset)
        self.show_result = show_result
        self.currently_remasking=False
        self.show_recalculated_fit=show_recalculated_fit

    def on_click(self,event):
        if (not event.inaxes is self.ax) or len(self.peaks)>2 or self.fig.canvas.toolbar.mode != "" or self.currently_remasking:
            return
        if len(self.new_points_plot)>2:
            try:
                self.new_points_plot[1].remove()
                self.new_points_plot.pop()
            except:
                pass
            try:
                self.new_points_plot[0].remove()
                self.new_points_plot.pop()
            except NotImplementedError:
                pass
            # self.fig.canvas.draw()
        c, v = int(event.xdata), self.y[int(event.xdata)]
        self.peaks.append((c, v))
        self.new_points_plot.append(self.ax.plot(c, v, 'x',color="green")[0])
        self.fig.canvas.draw()
        if len(self.peaks) == 2:
            self.currently_remasking=True
            if self.recalculated_fit:
                try:
                    self.recalculated_fit.remove()
                    self.recalculated_fit=None
                except NotImplementedError:
                    pass

            if self.peaks[0]>self.peaks[1]:
                temp=self.peaks[0]
                self.peaks[0]=self.peaks[1]
                self.peaks[1]=temp
                del temp

            self.results_of_remasking=self.func(self.peaks, self.y, self.x, self.params)
            masked_y, fitted_y, r2, t_R,sigma, D, R_h, t, p=self.results_of_remasking

            self.recalculated_fit=self.ax.plot(self.x, fitted_y,"--")[0]

            self.ax.set_ylim(min(*fitted_y,*self.y,*self.EIC_result.removed_dip_fitted), max(*fitted_y,*self.y,*self.EIC_result.removed_dip_fitted))
            self.abort_btn.setEnabled(True)
            self.continue_btn.setEnabled(True)
            self.reset_btn.setEnabled(False)
            self.abort_btn.clicked.connect(self.on_abort)
            self.continue_btn.clicked.connect(self.on_continue)

            self.fig.tight_layout()
            self.fig.canvas.draw_idle()
            self.show_recalculated_fit(*self.results_of_remasking)
            self.peaks=[]


    def on_abort(self):
        self.show_result()
        self.abort_btn.setEnabled(False)
        self.continue_btn.setEnabled(False)
        self.reset_btn.setEnabled(True)
        self.currently_remasking=False
        self.fig.canvas.draw_idle()

    def on_continue(self):
        masked_y, fitted_y, r2, t_R, sigma, D, R_h, t, p=self.results_of_remasking
        self.EIC_result.removed_dip = masked_y
        self.EIC_result.removed_dip_fitted = fitted_y
        self.EIC_result.r2 = r2
        self.EIC_result.tR = t_R
        self.EIC_result.sigma = sigma
        self.EIC_result.D = D
        self.EIC_result.Rh = R_h
        self.EIC_result.t = t
        self.EIC_result.p = p
        self.show_result()
        self.abort_btn.setEnabled(False)
        self.continue_btn.setEnabled(False)
        self.reset_btn.setEnabled(True)
        self.currently_remasking=False
        self.fig.canvas.draw_idle()

    def on_reset(self):
        self.EIC_result.removed_dip=self.EIC_result.original_removed_dip
        self.EIC_result.removed_dip_fitted=self.EIC_result.original_removed_dip_fitted
        self.EIC_result.r2=self.EIC_result.original_r2
        self.EIC_result.tR=self.EIC_result.original_tR
        self.EIC_result.sigma=self.EIC_result.original_sigma
        self.EIC_result.D=self.EIC_result.original_D
        self.EIC_result.Rh=self.EIC_result.original_Rh
        self.EIC_result.t=self.EIC_result.original_t
        self.EIC_result.p=self.EIC_result.original_p
        self.show_result()

    def clean_up(self):
        if self.cid:
            self.fig.canvas.mpl_disconnect(self.cid)
            self.cid=None
        try:
            self.reset_btn.clicked.disconnect(self.on_reset)
        except:
            pass

        try:
            self.abort_btn.clicked.disconnect(self.on_abort)
        except:
            pass

        try:
            self.continue_btn.clicked.disconnect(self.on_continue)
        except:
            pass

        if self.recalculated_fit:
            try:
                self.recalculated_fit.remove()
            except:
                pass

        for plot in self.new_points_plot:
            try:
                plot.remove()
            except:
                pass

        self.new_points_plot.clear()