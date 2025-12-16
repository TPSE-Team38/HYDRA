from typing import Callable
import matplotlib.pyplot as plt

class result():
    def __init__(self,y,x,params,fig,ax,func:Callable):
        self.peaks=[]
        self.y=y
        self.x=x
        self.changed=False
        self.new_points_plot=None
        self.params=params
        self.fig,self.ax=fig,ax
        self.cid=self.fig.canvas.mpl_connect("button_press_event",self.on_click)
        self.recalculated_fit = None
        self.func=func
        # self.recalculated_mask = None
        # self.recalculated_scatter = None

    def on_click(self,event):
        if not event.inaxes or len(self.peaks)>2:
            return
        c, v = int(event.xdata), event.ydata
        self.peaks.append((c, self.y[c]))
        self.new_points_plot=self.ax.plot(c, v, 'x')
        self.fig.canvas.draw()
        if len(self.peaks) == 2:
            if self.changed:
                self.recalculated_fit.remove()
                self.new_points_plot.remove()
                # self.recalculated_mask.remove()
                # self.recalculated_scatter.remove()
            if self.peaks[0]>self.peaks[1]:
                temp=self.peaks[0]
                self.peaks[0]=self.peaks[1]
                self.peaks[1]=temp
                del temp
            self.changed=True
            masked_y,fitted_y,r2,t_R,D,R_h=self.func(self.peaks, self.y, self.x, self.params)
            self.recalculated_fit=self.ax.plot(self.x, fitted_y,"--",label=f"recalculated_fit with r_2 score of {r2}\n t_R of {t_R}\n and diffusion coefficient of {D}\n and R_h of {R_h}")[0]
            # self.recalculated_mask=self.ax.plot(self.x, masked_y,"--",label="recalculated_mask")[0]
            # self.recalculated_scatter=self.ax.scatter(self.x, self.y,label="original")
            self.ax.set_ylim(min(*fitted_y,*self.y), max(*fitted_y,*self.y))
            self.fig.canvas.draw()
            self.ax.legend()
            self.peaks=[]
            plt.show()
