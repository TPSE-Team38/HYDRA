from typing import Callable

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib.artist import Artist
from matplotlib.lines import Line2D


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
    def __init__(self,y,x,params,fig,ax,func:Callable):
        self.peaks=[]
        self.y=y
        self.x=x
        self.new_points_plot=[]
        self.params=params
        self.fig,self.ax=fig,ax
        self.cid=self.fig.canvas.mpl_connect("button_press_event",self.on_click)
        self.recalculated_fit=None
        self.func=func
        # self.recalculated_mask = None
        # self.recalculated_scatter = None

    def on_click(self,event):
        if not event.inaxes or len(self.peaks)>2 or self.fig.canvas.toolbar.mode != "":
            return
        if len(self.new_points_plot)>1:
            self.new_points_plot[0].remove()
            self.new_points_plot[1].remove()
            self.new_points_plot.pop()
            self.new_points_plot.pop()
            self.fig.canvas.draw()
        c, v = int(event.xdata), self.y[int(event.xdata)]
        self.peaks.append((c, v))
        self.new_points_plot.append(self.ax.plot(c, v, 'x')[0])
        self.fig.canvas.draw()
        if len(self.peaks) == 2:
            if self.recalculated_fit:
                self.recalculated_fit.remove()
                self.recalculated_fit=None
                self.ax.legend()
                self.fig.canvas.draw()
                # self.recalculated_mask.remove()
                # self.recalculated_scatter.remove()

            if self.peaks[0]>self.peaks[1]:
                temp=self.peaks[0]
                self.peaks[0]=self.peaks[1]
                self.peaks[1]=temp
                del temp

            masked_y,fitted_y,r2,t_R,D,R_h,t,p=self.func(self.peaks, self.y, self.x, self.params)
            self.recalculated_fit=self.ax.plot(self.x, fitted_y,"--",label=f"recalculated_fit with r_2 score of {r2}\n t_R of {t_R}\n and diffusion coefficient of {D}\n and R_h of {R_h}")[0]
            # self.recalculated_mask=self.ax.plot(self.x, masked_y,"--",label="recalculated_mask")[0]
            # self.recalculated_scatter=self.ax.scatter(self.x, self.y,label="original")
            self.ax.set_ylim(min(*fitted_y,*self.y), max(*fitted_y,*self.y))
            # self.fig.canvas.draw()
            self.ax.legend()
            self.peaks=[]
            self.fig.canvas.draw_idle()