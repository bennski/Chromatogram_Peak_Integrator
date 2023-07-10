#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05-25-23

@author: bennski
"""

import plotly.offline as py
import plotly.graph_objects as go
import plotly.figure_factory as ff


import numpy as np
import pandas as pd
import scipy
import peakutils
import pylab as plb
import os

gates = [[11.58,12.35],[12.35,13.25],[13.45,13.87],[16.04,16.4],[16.73,17.15],[19.15,19.95],[19.95,22.15]] #uracil,thymine,cytosine,hypoxanthine,adenine,xanthine,guanine

noise_regions = [[11.83,12.35],[12.7,13.19],[13.75,13.87],[16.3,16.4],[17.,17.15],[19.5,19.95],[20.5,22.15]]

data = plb.loadtxt('sample_data2023.dat', skiprows=2)

indices = data[:,0]
time = data[:,1]
mag = data[:,2]

f = open('Nucleobase_peak_areas2023.dat', 'w')

f.write("Uracil\tUracil_errmin\tUracil_errmax\tThymine\tThymine_errmin\tThymine_errmax\tCytosine\tCytosine_errmin\tCytosine_errmax\tHypoxanthine\tHypoxanthine_errmin\tHypoxanthine_errmax\tAdenine\tAdenine_errmin\tAdenine_errmax\tXanthine\tXanthine_errmin\tXanthine_errmax\tGuanine\tGuanine_errmin\tGuanine_errmax\n")

for i in range(0,len(gates)):

    time_series = np.asarray(mag)

    area_under_first_peak = np.zeros(shape=(3))

    x = time
    time_series = time_series.tolist()

    diff_array = np.absolute(time - gates[i][0])
    gate_start_index = np.argmin((diff_array))
    if time[gate_start_index] < gates[i][0]:
        gate_start_index = gate_start_index + 1

    diff_array = np.absolute(time - gates[i][1])
    gate_end_index = np.argmin((diff_array))
    if time[gate_end_index] > gates[i][1]:
        gate_end_index = gate_end_index - 1

    diff_array = np.absolute(time - noise_regions[i][0])
    noise_start_index = np.argmin((diff_array))
    if time[noise_start_index] < noise_regions[i][0]:
        noise_start_index = noise_start_index + 1

    diff_array = np.absolute(time - noise_regions[i][1])
    noise_end_index = np.argmin((diff_array))
    if time[noise_end_index] > noise_regions[i][1]:
        noise_end_index = noise_end_index - 1

    noise_floor_lower = min(mag[gate_start_index:gate_end_index+1])
    if noise_floor_lower < 1:
        noise_floor_lower = 1

    peak_index = np.argmax(mag[gate_start_index:gate_end_index+1])


    noise_average = np.average(mag[noise_start_index:noise_end_index+1])
    noise_std = np.std(mag[noise_start_index:noise_end_index+1])
    noise_plus_std = noise_average + noise_std
    noise_minus_std = noise_average - noise_std
    print('Peak height: %s' % (np.round(np.max(mag[gate_start_index:gate_end_index+1]),2)))
    print('Noise avg: %s +- %s' % (round(noise_average,2), round(noise_std,2)))
    SNR_Val = round(np.max(mag[gate_start_index:gate_end_index+1])/noise_average,1)
    print('SNR: %s' % (SNR_Val))
    if noise_average < 0:
        noise_average = 0
    if noise_plus_std < 0:
        noise_plus_std = 0
    if noise_minus_std < 0:
        noise_minus_std = 0


    #NEW_AVERAGE_CALC
    baseline_values = np.ones(shape=(len(time_series)))*noise_average
    baseline_values = baseline_values.tolist()

    start_index = gate_start_index
    end_index = gate_end_index+1
    for j in range(gate_start_index+peak_index,gate_start_index,-1):
        if j < gate_start_index+peak_index:
            if mag[j] < noise_average or mag[j] > mag[j+1]:
                start_index = j
                break
    for j in range(gate_start_index+peak_index,gate_end_index+1):
        if j > gate_start_index+peak_index:
            if mag[j] < noise_average or mag[j] > mag[j-1]:
                end_index = j
                break

    x_start = np.interp(noise_average,time_series[start_index:start_index+2],x[start_index:start_index+2])

    x_finish = np.interp(noise_average,time_series[end_index-1:end_index+1][::-1],x[end_index-1:end_index+1][::-1])

    x_vals = np.append(x_start,x[start_index+1:end_index])
    x_vals = np.append(x_vals,x_finish)

    y_vals = np.append(noise_average,time_series[start_index+1:end_index])
    y_vals = np.append(y_vals,noise_average)

    area_under_first_peak[0] = np.trapz(y_vals, x_vals) - np.trapz(baseline_values[start_index:end_index+1], x_vals)
    print('Average peak integration %s' % (np.round(area_under_first_peak[0],2)))

    rev_baseline_values = baseline_values[start_index-1:end_index+1]
    rev_baseline_values.reverse()
    area_x = np.append(x_vals,x_vals[::-1])
    area_y = np.append(y_vals,rev_baseline_values)

    #Uncertainty_CALC_min
    baseline_values_min = np.ones(shape=(len(time_series)))*(noise_plus_std)
    baseline_values_min = baseline_values_min.tolist()

    start_index = gate_start_index
    end_index = gate_end_index+1
    for j in range(gate_start_index+peak_index,gate_start_index,-1):
        if j < gate_start_index+peak_index:
            if mag[j] < (noise_plus_std) or mag[j] > mag[j+1]:
                start_index = j
                break
    for j in range(gate_start_index+peak_index,gate_end_index+1):
        if j > gate_start_index+peak_index:
            if mag[j] < (noise_plus_std) or mag[j] > mag[j-1]:
                end_index = j
                break


    x_start = np.interp((noise_plus_std),time_series[start_index:start_index+2],x[start_index:start_index+2])

    x_finish = np.interp(((noise_plus_std)),time_series[end_index-1:end_index+1][::-1],x[end_index-1:end_index+1][::-1])

    x_vals = np.append(x_start,x[start_index+1:end_index])
    x_vals = np.append(x_vals,x_finish)

    y_vals = np.append((noise_plus_std),time_series[start_index+1:end_index])
    y_vals = np.append(y_vals,(noise_plus_std))

    area_under_first_peak[1] = np.trapz(y_vals, x_vals) - np.trapz(baseline_values_min[start_index:end_index+1], x_vals)
    print('Min peak integration %s' % (np.round(area_under_first_peak[1],2)))


    #Uncertainty_CALC_max
    baseline_values_max = np.ones(shape=(len(time_series)))*(noise_minus_std)
    baseline_values_max = baseline_values_max.tolist()

    start_index = gate_start_index
    end_index = gate_end_index+1
    for j in range(gate_start_index+peak_index,gate_start_index,-1):
        if j < gate_start_index+peak_index:
            if mag[j] < (noise_minus_std) or mag[j] > mag[j+1]:
                start_index = j
                break
    for j in range(gate_start_index+peak_index,gate_end_index+1):
        if j > gate_start_index+peak_index:
            if mag[j] < (noise_minus_std) or mag[j] > mag[j-1]:
                end_index = j
                break

    x_start = np.interp((noise_minus_std),time_series[start_index:start_index+2],x[start_index:start_index+2])

    x_finish = np.interp((noise_minus_std),time_series[end_index-1:end_index+1][::-1],x[end_index-1:end_index+1][::-1])

    x_vals = np.append(x_start,x[start_index+1:end_index])
    x_vals = np.append(x_vals,x_finish)

    y_vals = np.append((noise_minus_std),time_series[start_index+1:end_index])
    y_vals = np.append(y_vals,(noise_minus_std))

    area_under_first_peak[2] = np.trapz(y_vals, x_vals) - np.trapz(baseline_values_max[start_index:end_index+1], x_vals)
    print('Max peak integration %s' % (np.round(area_under_first_peak[2],2)))


    if i == 0:
        molecule = "Uracil"
    if i == 1:
        molecule = "Thymine"
    if i == 2:
        molecule = "Cytosine"
    if i == 3:
        molecule = "Hypoxanthine"
    if i == 4:
        molecule = "Adenine"
    if i == 5:
        molecule = "Xanthine"
    if i == 6:
        molecule = "Guanine"


    print('The average area for %s is %s - %s + %s' % (molecule,np.round(area_under_first_peak[0],2),np.round((area_under_first_peak[0]-area_under_first_peak[1]),2),np.round((area_under_first_peak[2]-area_under_first_peak[0]),2)))
    f.write(str(np.round(np.average(area_under_first_peak[0]),2)))
    f.write("\t")
    f.write(str(np.round((area_under_first_peak[0]-area_under_first_peak[1]),2)))
    f.write("\t")
    f.write(str(np.round((area_under_first_peak[2]-area_under_first_peak[0]),2)))
    f.write("\t")


    trace = go.Scatter(
        x=x,
        y=time_series,
        mode='lines',
        marker=dict(
            color='#B292EA',
            size=5,
        ),
        line=dict(width=5,
        ),
        name='Signal'
    )

    trace2 = go.Scatter(
        x=x,
        y=baseline_values,
        mode='lines',
        marker=dict(
            size=3,
            color='#444444',
        ),
        name='Average Noise Baseline'
    )

    trace2min = go.Scatter(
        x=x,
        y=baseline_values_min,
        mode='lines',
        marker=dict(
            size=3,
            color='#bbbbbb',
            symbol="diamond",
        ),
        name='Baseline + Sigma'
    )

    trace2max = go.Scatter(
        x=x,
        y=baseline_values_max,
        mode='lines',
        marker=dict(
            size=3,
            color='#bbbbbb',
            symbol="square",
        ),
        name='Baseline - Sigma'
    )

    trace3 = go.Scatter(
        x=area_x,
        y=area_y,
        mode='lines+markers',
        marker=dict(
            size=14,
            color='rgb(255,0,0)',
        ),
        name='Peak Outline'
    )

    trace4 = go.Scatter(
        x=np.ones(shape=(1000))*gates[i][0],
        y=np.linspace(-20000,20000,1000),
        mode='lines',
        marker=dict(
            size=3,
            color='#4334eb',
        ),
        name='Gate Start'
    )

    trace5 = go.Scatter(
        x=np.ones(shape=(1000))*gates[i][1],
        y=np.linspace(-20000,20000,1000),
        mode='lines',
        marker=dict(
            size=3,
            color='#4334eb',
        ),
        name='Gate End'
    )

    trace6 = go.Scatter(
        x=np.linspace(noise_regions[i][0],noise_regions[i][1],len(time_series[noise_start_index:noise_end_index+1])),
        y=time_series[noise_start_index:noise_end_index+1],
        mode='markers',
        marker=dict(
            size=8,
            color='#eded18',
            symbol="diamond"
        ),
        name='Noise Region'
    )


    layout = go.Layout()

    trace_data = [trace, trace2, trace2min, trace2max, trace3,trace4,trace5,trace6]
    fig = go.Figure(data=trace_data, layout=layout)

    fig.update_layout(
    hoverlabel=dict(
        font_size=50,
        font_family="Rockwell"
    ),
    yaxis = dict(
        tickfont = dict(size=40)),
    xaxis = dict(
        tickfont = dict(size=40)),
    legend=dict(title_font_family="Times New Roman",
                              font=dict(size= 28))
)

    if i==2:
        fig.add_annotation(x=13.75, y=800,
            text="SNR = "+str(SNR_Val),
            showarrow=False,
            yshift=10,
            font_size=30)

    #fig.update_yaxes(type="log")
    py.iplot(fig, filename='Nucleobase-peak-integration')

f.close()