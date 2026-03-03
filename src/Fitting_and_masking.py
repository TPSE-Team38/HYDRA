import numpy as np
from scipy.optimize import curve_fit
from .Calculations import gaus,new_gauss_from_jonathan
import scipy.signal as sig
from scipy import interpolate
'''might be needed'''

'''masking'''
##-----------------------------##
#finds the peaks of the smoothed graph
#finds the smallest value between first and last peak: min
#then sets the mid point of interpolation to double of: biggest peak - min
#interpolates between first peak, mid point, last peak quadratically
def dip_detect_correct(y):
    """
    finds the peaks of the smoothed graph
    finds the smallest value between first and last peak: min
    then sets the mid point of interpolation to double of: biggest peak - min
    interpolates between first peak, mid point, last peak quadratically
    y: smoothed values of the y-axis
    Returns graph with dip replaced with interpolated values"""

    #get middle point by adding the highest peak to the difference between the highest peek and the lowest point
    peak, _ = sig.find_peaks(y, height=np.mean(y))
    max_peak=max(y[peak])
    minPoint = max_peak + (max_peak - min(y[peak[0]:peak[-1]]))

    pointsX=np.array([peak[0] - 1,
                      peak[-1] - ((peak[-1] - peak[0]) / 2),
                      peak[-1] + 1])
    pointsY=np.array([y[peak[0] - 1]
                         , minPoint
                         ,y[peak[-1] + 1]])
    #interpolate using first and last peak in x-axis and middle-point
    predict = interpolate.interp1d(pointsX,pointsY,kind='quadratic')

    #connect the interpolated values to the tails on either side of the left and right peak
    interpolated = np.concatenate((y[:(peak[0] - 1)], predict(np.arange(1,len(y))[peak[0] - 1:peak[len(peak) - 1] + 1]),
                                           y[peak[len(peak) - 1] + 1:]))

    return interpolated

def mask_region(y,left,right):
    return np.concatenate((y[:left],[np.nan]*(right-left),y[right:]))

#doesn't work
'''
def new_masking(y:np.ndarray[float],x:np.ndarray[float],params:np.ndarray[float]):
    points=set([x for x in range(1,len(y)-1) if abs(x-np.argmax(y)) < 3 ])
    rng=np.random.default_rng()
    highest_r2=-1
    highest_r2_point=(-1,-1)
    highest_midpoint=-1
    besttau=-1
    leftDone=False
    rightDone=False
    while len(points)>0:
        if besttau>=1.25:
            break
        mid = rng.choice(list(points))
        leftP, rightP = mid - 1, mid + 1
        points.remove(mid)
        if abs(mid-np.argmax(y))>=1:
            print(f"skipped point {mid}")
            continue
        print(f"testing point {mid}")
        saveL, saveR = leftP, rightP
        leftDone,rightDone=False,False
        while not leftDone or not rightDone:
            # print(i,highest_r2_point)
            if not leftDone:
                if rightP>=len(y):
                    if leftP<=0 or abs(mid-leftP)>len(y)/12:
                        leftDone=True
                    else:
                        leftP-=1
                    rightP=saveR
                else:
                    rightP+=1
            else:
                if not rightDone:
                    if leftP>0:
                        if rightP>=len(y) or abs(mid-rightP)>len(y)/12:
                            rightDone=True
                        else:
                            rightP+=1
                        leftP=saveL
                    else:
                        leftP+=1
            print("fu")
            masked_y,fitted,r2,t_R,D,R_h=recalculate([(leftP,0),(rightP,0)],y,x,params)
            # print(params[0])
            t=tau(params[0],params[3],params[1],params[4],R_h)
            # masked_y=mask_region(y,leftP,rightP)
            # mask=~np.isnan(masked_y)
            fitted_y,sigma=different_approach_gaus(masked_y,x,len(y)/2)
            # r_2=r2_score(masked_y[mask],fitted_y[mask])
            if t>besttau:
                besttau=t
                # highest_r2=r_2
                highest_r2_point=(leftP,rightP)
                highest_midpoint=mid
                print(besttau,",",highest_midpoint)
            if besttau >= 1.25:
                break

        print(f"done testing point {mid}")

    return mask_region(y,*highest_r2_point),highest_midpoint
'''
def get_peaks(y):
    #try getting peaks beginning from left then beginning from right (sigmoidal fit)

    maxDiff=0
    for i in range(0,len(y)-2):
        maxDtemp=max(abs(y[i]-y[i+1]),abs(y[i]-y[i+2]))
        if maxDtemp>maxDiff or maxDtemp>maxDiff:
            maxDiff=maxDtemp
    '''
    maxDiff=0
    for i in range(0,len(y)-1):
        maxDtemp=abs(y[i]-y[i+1])
        if maxDtemp>maxDiff or maxDtemp>maxDiff:
            maxDiff=maxDtemp
    '''
    leftP=0
    for i,x in enumerate(y):
        if x>y[leftP]:
            leftP=i
        elif abs(x-y[leftP])>maxDiff:
            break
    rightP=0
    for i in reversed(range(1,len(y))):
        if y[i]>y[rightP]:
            rightP=i
        elif abs(y[i]-y[rightP])>maxDiff:
            break
    # y_new=np.concatenate((y[:leftP],[np.nan]*len(y[leftP:rightP]),y[rightP:]))
    actual_right = -1
    actual_left = -1
    right_side = y[rightP:]
    left_side = y[:leftP+1]
    done = False
    boundL, boundR = 0.94, 0.95
    while not done:
        for (i, x) in enumerate(right_side):
            if boundL * y[rightP] <= x <= boundR * y[rightP] :
                actual_right = rightP+i
                break
        if actual_right == -1:
            boundL -= 0.0001
        else:
            done = True
    done = False
    boundL, boundR = 0.94, 0.95
    while not done:
        for (i, x) in enumerate(left_side):
            if boundL * left_side[leftP] <= x <= boundR * left_side[leftP]:
                actual_left = i
                break
        if actual_left == -1:
            boundL -= 0.0001
        else:
            done = True
    y_new=np.concatenate((y[:actual_left+1],[np.nan]*len(y[actual_left+1:actual_right]),y[actual_right:]))
    return y_new,leftP+((rightP-leftP)/2)

#not being used
def get_2_peaks(y:np.ndarray[float],y_nonsmoothed)->np.ndarray[float]:
    """
    Returns y-values with the found noise created by ion suppression removed(replaced by np.nan)
    """

    y_smoothed=sig.savgol_filter(y_nonsmoothed,window_length=11,polyorder=2,mode='mirror')
    # plt.figure()
    # plt.plot(np.arange(1,len(y_smoothed)+1),y_smoothed)
    # plt.savefig("fig.png")
    peaks,_=sig.find_peaks(y_smoothed,height=y_smoothed.mean())
    x_c=peaks[0]+np.argmin(y[peaks[0]:peaks[-1]])

    right_peak=np.argmax(y_nonsmoothed[x_c+1:])
    left_peak=np.argmax(y_nonsmoothed[:x_c])
    actual_right=-1
    actual_left=-1
    right_side=y_nonsmoothed[right_peak+x_c:]
    left_side=y_nonsmoothed[:left_peak]
    done=False
    boundL,boundR=0.94,0.96
    while not done:
        for (i,x) in enumerate(right_side):
            if boundL*right_side[right_peak]<= x <=boundR*right_side[right_peak]:
                actual_right=x_c+i
                break
        if actual_right==-1:
            boundL-=0.01
        else:
            done=True
    done=False
    boundL,boundR=0.94,0.96
    while not done:
        for (i,x) in enumerate(left_side):
            if boundL*left_side[left_peak-1]<= x <=boundR*left_side[left_peak-1]:
                actual_left=i
                break
        if actual_left==-1:
            boundL-=0.01
        else:
            done=True
    maxPeak=max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right])
    if (maxPeak-y_nonsmoothed[x_c])/maxPeak>0.3:
        new_y=np.concatenate((y_nonsmoothed[:actual_left],[np.nan]*len(y[actual_left:actual_right]),y_nonsmoothed[actual_right:]))
        # print((max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right])-y_nonsmoothed[x_c])/max(y_nonsmoothed[actual_left],y_nonsmoothed[actual_right]))
    else:
        new_y=y_nonsmoothed.copy()
        fit=gaussian_fit(new_y,np.arange(1,len(new_y)+1))
        inter=np.intersect1d(y,fit)
        # print(inter)

    return new_y,x_c

################################################################
'''fitting'''
##-----------------------------##
#different_approach_gaus_jonathan
def different_approach_gaus_jonathan(y:np.ndarray[float],x:np.ndarray[float],xc):
    """
    fits the curve of the tails created by ion suppression removal
    """
    mask=~np.isnan(y)
    x_fit=x[mask]
    y_fit=y[mask]
    # p0 = [min(y_fit), xc,2.355*(np.std(y_fit)), (y_fit.max() - y_fit.min())]
    params,_=curve_fit(new_gauss_from_jonathan,x_fit,y_fit,maxfev=10000)
    return new_gauss_from_jonathan(x,*params),(params[2]/2)

def gaussian_fit(y:np.ndarray[float],x:np.ndarray[float],xc):
    """
    fits the curve of the tails created by ion suppression removal
    """
    mask=~np.isnan(y)
    x_fit=x[mask]
    y_fit=y[mask]
    mu=np.mean(y_fit)
    sum_squared_deviation=sum([(p-mu)**2 for p in y_fit])
    variance=sum_squared_deviation/len(y_fit)-1
    sigma=np.sqrt(variance)
    # p0=[max(y)-np.median(y_fit,axis=0),xc,sigma,np.median(y_fit,axis=0)]
    try:
        params,_=curve_fit(gaus,x_fit,y_fit,maxfev=10000,method='lm')
    except:
        params=[0,0,0,0]
    return gaus(x,*params),params[2]

