# functions for helping analysis ICP-MS data

# imports
import numpy as np
import re
from scipy.optimize import curve_fit
from math import nan
from scipy import interpolate
import matplotlib.pyplot as plt
import pickle
import csv
import math

#----------------Read ICP-MS File--------------------------------------------#

# extracts ppb data (mean, sd) from csv file
#WARNING: I DO NOT YET KNOW HOW SUBJECT TO CHANGE THE FILE STRUCTURE IS
def read_icpms_file(file_name, has_rsd = True):
    with open(file_name, newline='') as csvfile:
        rows= csv.reader(csvfile, delimiter=',', quotechar='"')
        els = rows.__next__()
        el_start = 9 # HARDCODED FOR NOW; COULD VARY FROM RUN TO RUN
        samp_dic = {}
        for row in rows:
            if row[6] == '':
                el_dic = {}
                samp_name = row[7]
                step = 2
                if not has_rsd:
                    step = 1
                for i in range(el_start,len(row),step):
                    try:
                        mean_val = float(row[i])
                        sd_val = 0
                        if has_rsd:
                            sd_val = float(row[i+1])
                        el_dic[els[i]] = (mean_val, sd_val)
                    except Exception:
                        try:
                            mean_val = float(row[i])
                            el_dic[els[i]] = (mean_val, 0)
                        except Exception:
                            continue
                samp_dic[samp_name] = el_dic
    return samp_dic

def read_icpms_file_in_order(file_name, has_rsd = True):
    with open(file_name, newline='') as csvfile:
        rows= csv.reader(csvfile, delimiter=',', quotechar='"')
        els = rows.__next__()
        el_start = 9 # HARDCODED FOR NOW; COULD VARY FROM RUN TO RUN
        samp_names = []
        samp_data = []
        for row in rows:
            if row[6] == '':
                el_dic = {}
                samp_name = row[7]
                step = 2
                if not has_rsd:
                    step = 1
                for i in range(el_start,len(row),step):
                    try:
                        mean_val = float(row[i])
                        sd_val = 0
                        if has_rsd:
                            sd_val = float(row[i+1])
                        el_dic[els[i]] = (mean_val, sd_val)
                    except Exception:
                        try:
                            mean_val = float(row[i])
                            el_dic[els[i]] = (mean_val, 0)
                        except Exception:
                            continue
                samp_names.append(samp_name)
                samp_data.append(el_dic)
    return samp_names, samp_data

def read_cal_curve(file_name, has_rsd = True):
    with open(file_name, newline='') as csvfile:
        rows= csv.reader(csvfile, delimiter=',', quotechar='"')
        els = rows.__next__()
        el_start = 9 # HARDCODED FOR NOW; COULD VARY FROM RUN TO RUN
        samp_dic = {}
        for row in rows:
            if row[5] == 'CalStd':
                el_dic = {}
                samp_name = row[7]
                step = 2
                if not has_rsd:
                    step = 1
                for i in range(el_start,len(row),step):
                    try:
                        mean_val = float(row[i])
                        sd_val = 0
                        if has_rsd:
                            sd_val = float(row[i+1])
                        el_dic[els[i]] = (mean_val, sd_val)
                    except Exception:
                        try:
                            mean_val = float(row[i])
                            el_dic[els[i]] = (mean_val, 0)
                        except Exception:
                            continue
                samp_dic[samp_name] = el_dic
    return samp_dic


# extracts uM measurements from dilution info and ppb data
def get_uM_meas(samp, ree_info, dil):
    m_mw = re.search(r'(\d+)', ree_info)
    mw = float(m_mw.group(1))
    mean_val = samp[ree_info][0] / dil / mw
    sd_val = samp[ree_info][1] / 100
    return (mean_val,sd_val)

def max_var_analysis(pts1, pts2):
    # calculates center points
    mean1 = np.mean(pts1)
    mean2 = np.mean(pts2)
    
    # shifts points such that the means are zero
    pts1_adj = np.array([pts1[i] - mean1 for i in range(len(pts1))])
    pts2_adj = np.array([pts2[i] - mean2 for i in range(len(pts2))])
    
    # finds how much we need to rotate axes to maximize variation along one axis
    theta = 0.5 * math.atan(np.sum(2 * pts1_adj * pts2_adj) / np.sum(pts1_adj**2 - pts2_adj**2))
    
    # rotates shifted points by theta
    pts1_rot = pts1_adj * math.cos(theta) + pts2_adj * math.sin(theta)
    pts2_rot = pts2_adj * math.cos(theta) - pts1_adj * math.sin(theta)
    
    # calculates standard deviation on each new axis
    axis1_sd = np.std(pts1_rot)
    axis2_sd = np.std(pts2_rot)
    
    # finds points needed for errorbars
    err1_xs = [-axis1_sd * math.cos(theta) + mean1, axis1_sd * math.cos(theta) + mean1]
    err1_ys = [axis1_sd * math.sin(-theta) + mean2, -axis1_sd * math.sin(-theta) + mean2]
    err2_xs = [-axis2_sd * math.sin(-theta) + mean1, axis2_sd * math.sin(-theta) + mean1]
    err2_ys = [-axis2_sd * math.cos(theta) + mean2, axis2_sd * math.cos(theta) + mean2]
    
    return err1_xs, err1_ys, err2_xs, err2_ys

def fit_avg_curve(pts1, pts2, num_steps_per_pt = 2, min_pts_per_step = 2):
    # sorts points by pt1
    idxes = np.argsort(pts1)
    pts1 = np.array([pts1[i] for i in idxes])
    pts2 = np.array([pts2[i] for i in idxes])
    
    # calculates step interval and calculates pt2 averages at each step
    steps, step = np.linspace(pts1[0], pts1[-1], len(pts1) // num_steps_per_pt, True, True)
    avgs = np.zeros(len(steps)) # average of pt2 values
    num_cont = np.zeros(len(steps)) # number of points contributing to the averages
    for i in range(len(steps)):
        idxes1 = np.argwhere(np.abs(pts1 - steps[i]) <= step)
        avgs[i] = np.median([pts2[j] for j in idxes1])
        num_cont[i] = len(idxes1)
    
    # trims to only include areas with sufficient data
    long_arr = []
    cur_arr = []
    for i in range(len(steps)):
        if num_cont[i] >= min_pts_per_step:
            cur_arr.append(i)
        else:
            if len(cur_arr) > len(long_arr):
                long_arr = cur_arr
            cur_arr = []
    if len(cur_arr) > len(long_arr):
        long_arr = cur_arr
        
    return steps[long_arr[0]:(long_arr[-1] + 1)], avgs[long_arr[0]:(long_arr[-1] + 1)]

def fit_line_curve(pts1, pts2, lim = 0.1, deg = 1):
    idxes = np.argsort(pts1)
    start = int(len(pts1) * lim)
    end = int(len(pts1) * (1 - lim))
    pts1 = np.array([pts1[i] for i in idxes[start:end]])
    pts2 = np.array([pts2[i] for i in idxes[start:end]])
    line_fit = np.polyfit(pts1, pts2, deg)
    return line_fit, pts1[0], pts1[-1] # returns best fit, range of fit