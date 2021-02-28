#-------------------------------------------------------------------------------
# Name:        pythreshold_support.py
# Purpose:     Functions support for processin of pythresholdsERA5_yemen.py
#-------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shapefile as shp
from datetime import datetime, timedelta

def as_stride(arr, sub_shp, stride):
    """Get a strided sub-matrices view"""
    s0, s1 = arr.strides[:2]
    m1, n1 = arr.shape[:2]
    m2, n2 = sub_shp
    view_shp = (1+(m1-m2)//stride[0], 1+(n1-n2)//stride[1], m2, n2) + arr.shape[2:]
    strides = (stride[0]*s0, stride[1]*s1, s0, s1) + arr.strides[2:]
    subs = np.lib.stride_tricks.as_strided(arr, view_shp, strides=strides)
    return subs

def pooling_overlap(mat, ksize, stride=None, method='max', pad=False):
    """Overlaping pooling"""
    m, n = mat.shape[:2]
    ky, kx = ksize
    if stride is None:
        stride = (ky, kx)
    sx, sy = stride

    _ceil = lambda x, y : int(np.ceil(x/float(y)))

    if pad:
        ny = _ceil(m, sy)
        nx = _ceil(n, sx)
        size = ((ny-1)*sy+ky, (nx-1)*sx+kx) + mat.shape[2:]
        mat_pad = np.full(size, np.nan)
        mat_pad[:m, :n, ...] = mat
    else:
        mat_pad = mat[:(m-ky)//sy*sy+ky, :(n-kx)//sx*sx+kx, ...]

    view = as_stride(mat_pad, ksize, stride)

    if method == 'max':
        result = np.nanmax(view, axis=(2,3))
    elif method == 'min':
        result = np.nanmin(view, axis=(2,3))
    elif method == 'mean':
        result = np.nanmean(view, axis=(2,3))
    else:
        result = None
        raise Exception("Function not implemented.")

    return result

def pooling(mat, ksize, method='max'):
    """Pooling with single stride."""
    return pooling_overlap(mat, ksize, stride=(1,1), method=method)

def thres_step(d_arr, d_frm, threshold, method='max'):
    """Maximum value in a centered 7X7 box exceeds a threshold value."""
    counter = np.zeros(len(d_frm), dtype=np.int0)
    ksize = (7,7)
    val = ksize[0] // 2
    dim = len(d_arr.shape)
    assert dim in (3,4), "DataArray must be 3D or 4D."
    mom = d_arr.values if dim==3 else np.moveaxis(d_arr.values, 1, -1)
    dspoll = np.array([pooling(mom[i,...], ksize, method=method) for i in range(mom.shape[0])])
    for i in range(len(d_frm)):
        a = d_frm["LON_index"][i]
        b = d_frm["LAT_index"][i]
        c = d_frm["TIME_index"][i]
        v = dspoll[c, b-val, a-val] if dim==3 else dspoll[c, b-val, a-val, ...].mean()
        if v > threshold:
            counter[i] = 1
        #print(d_frm.index[i])
    d_frm['sub'] = counter
    return d_frm

def center_step(d_arr, d_frm):
    """Minimum value at center in a centered 7X7 box."""
    counter = np.zeros(len(d_frm), dtype=np.int0)
    ksize = (7,7)
    val = ksize[0] // 2
    dim = len(d_arr.shape)
    assert dim in (3,), "DataArray must be 3D."
    mom = d_arr.values
    dspoll = np.array([pooling(mom[i, ...], ksize, 'min') for i in range(mom.shape[0])])
    for i in range(len(d_frm)):
        a = d_frm["LON_index"][i]
        b = d_frm["LAT_index"][i]
        c = d_frm["TIME_index"][i]
        cond = np.argwhere(mom[c,b-val:b+val+1,a-val:a+val+1] == dspoll[c,b-3,a-3])
        if np.all(cond == np.array([val,val]), axis=1).any():
            counter[i] = 1
        #print(d_frm.index[i])
    d_frm['sub'] = counter
    return d_frm

def levels_positive_step(d_arr, d_frm):
    counter = np.zeros(len(d_frm), dtype=np.int0)
    ksize = (7,7)
    val = ksize[0] // 2
    dim = len(d_arr.shape)
    assert dim in (4,), "DataArray must be 4D."
    mom = np.moveaxis(d_arr.values, 1, -1)
    dspoll = np.array([pooling(mom[i,...], ksize, 'mean') for i in range(mom.shape[0])])
    for i in range(len(d_frm)):
        a = d_frm["LON_index"][i]
        b = d_frm["LAT_index"][i]
        c = d_frm["TIME_index"][i]
        if np.alltrue(dspoll[c, b-val, a-val, :]>0):
            counter[i] = 1
        #print(d_frm.index[i])
    d_frm['sub'] = counter
    return d_frm

def levels_compare_step(da_left, da_right,  d_frm):
    counter = np.zeros(len(d_frm), dtype=np.int0)
    ksize = (7,7)
    val = ksize[0] // 2
    dim = len(da_left.shape)
    assert dim in (3,4), "DataArray must be 3D or 4D."
    mom = da_left.values if dim==3 else da_left.values[:,0,...]
    dspoll_left = np.array([pooling(mom[i,...], ksize, 'mean') for i in range(mom.shape[0])])
    dim = len(da_right.shape)
    assert dim in (3,4), "DataArray must be 3D or 4D."
    mom = da_right.values if dim==3 else da_right.values[:,3,...]
    dspoll_right = np.array([pooling(mom[i,...], ksize, 'mean') for i in range(mom.shape[0])])
    for i in range(len(d_frm)):
        a = d_frm["LON_index"][i]
        b = d_frm["LAT_index"][i]
        c = d_frm["TIME_index"][i]
        if dspoll_left[c, b-3, a-3] > dspoll_right[c, b-3, a-3]:
            counter[i] = 1
        #print(d_frm.index[i])
    d_frm['sub'] = counter
    return d_frm

def fill_results(d_frm, d_arr1, d_arr2, d_arr3):
    fillerdtm = np.zeros((1, len(d_frm)), dtype='datetime64[ns]')
    fillerint = np.zeros((1, len(d_frm)), dtype=np.int0)
    fillerdbl = np.zeros((2, len(d_frm)), dtype=np.float64)
    for i in range(len(d_frm)):
        a = d_frm["LON_index"][i]
        b = d_frm["LAT_index"][i]
        c = d_frm["TIME_index"][i]
        fillerdtm[0, i] = d_arr1.values[c]
        fillerdbl[0, i] = d_arr2.values[c, b, a]
        fillerdbl[1, i] = d_arr3.values[c, b, a]
        fillerint[0, i] = i
    d_frm['TIME'] = np.int0((fillerdtm[0, :]-np.datetime64('1900-01-01T00:00:00Z')) / np.timedelta64(1, 'h'))
    d_frm['Eye_Press'] = fillerdbl[0, :]
    d_frm['Eye_Press2'] = fillerdbl[1, :]
    d_frm['ID'] = fillerint[0, :]
    d_frm['TIME1'] = fillerdtm[0, :]
    del fillerdtm
    del fillerint
    del fillerdbl
    return d_frm

def plot_bars_filters(df):
    ax = df.plot(x='Step', y='No.P', kind='bar', legend=None, width=0.8, color='grey')
    ax.set_title('No. of guessed cyclones eye after each filter')
    ax.set_ylabel('No.P')
    ax.xaxis.set_tick_params(rotation=0)
    rects = ax.patches
    for rec, val in zip(rects, df['No.P']):
        h = rec.get_height()
        ax.text(rec.get_x()+rec.get_width()/2, h+5, str(val), ha='center', va='bottom')
    #ax.grid(b=True, color='grey', linestyle='-.', linewidth=0, alpha=0.2)
    ax.grid(color='grey', linewidth=0.5, linestyle='-.')
    ax.set_axisbelow(True)
    plt.pause(1)

def plot_poins_pressure(df):
    ax = df.plot(x='TIME', y='Eye_Press', kind='scatter', legend=None, color='black')
    ax.set_xticklabels([datetime.strftime(date, '%d %b %Y\n%H:%M:%S') for date in df['TIME1']])
    ax.xaxis.set_tick_params(rotation=45)
    ax.grid(color='grey', linewidth=0.5, linestyle='-.')
    ax.set_axisbelow(True)
    plt.pause(1)

def plot_shape(filename):
    shapefile = shp.Reader(filename)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for shape in shapefile.shapes():
        points = np.array(shape.points)
        intervals = list(shape.parts)+[len(shape.points)]
        ymin, xmin = np.min(points, axis=0)
        ymax, xmax = np.max(points, axis=0)
        for (i, j) in zip(intervals[:-1], intervals[1:]):
            pts = points[i:j]
            if shapefile.shapeType == shp.POLYGON:
                ax.plot(*zip(*np.append(pts, pts[:1], axis=0)),'r')
            elif shapefile.shapeType == shp.POLYLINE:
                ax.plot(*zip(*pts),'b')
            elif shapefile.shapeType == shp.POINT:
                pass

    limits = ax.axis('off')
    ax.set_aspect(1)
    fig.tight_layout(pad=0.5)
    plt.pause(1)


def main():
    pass

if __name__ == '__main__':
    main()
