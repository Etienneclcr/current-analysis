# -*- coding: utf-8 -*-
"""
Created on Fri May 21 13:39:34 2021

@author: Etienne.CLAUSS
"""


import pandas as pd, winsound, numpy as np
from ToolKit_PSC import TraceManager_PSC, Analysor_PSC_2_0, PlotFactory_PSC_2_0


folder_exp = r"C:\Etienne.CLAUSS\Lab\Exp\ephy\BA\EPSC\Retrobeads_CeM_EPSC_ren"


'''                   Define parameters of the analysis                   '''
        
#Define parameters of the analysis & create empty lists
drug, event_type, phenotype = (folder_exp.split('\\')[-1], 
                               folder_exp.split('\\')[-2], 
                               folder_exp.split('\\')[-3])
drug_arrival, exp_duration, pol, bin_10s, bin_1s = 270, 1200, -1, 10, 1
chosen_wind, wind_bl, wind_dr, wind_wsh, wind_length = True, 120, 400, 900, 180
Paired = True
colorz = ('blue', 'limegreen', 'red', 'darkorange', 'lawngreen')

#Create folders to stock the analysis
FM = TraceManager_PSC.FileManager(folder_exp)
(folder_exp_analysis, individual_traces, 
 stock , writer) = FM.FileMaker(wind_length, drug)

#Find the files containing the data
folder_exp_list, n_cell = FM.FileSeeker()

#Create the Output time course dictionnary 
Stock = {x:[] for x in ('Cell ID', 'Binary', 'Spike Total', 'Bin 10s', 'Bin 1s', 
                        'bl_spikes', 'Ratio 10s', 'Indexes', 'Ecdf', 'Ampli')}

#Create a loop to analyse the cells one by one
for file_name in folder_exp_list:
    file = f'{folder_exp}\{file_name}'
    print('='*30, '\n'*2, file_name, '\n'*2, '='*30)
    
    '''                          Pre-Analysis part                     '''
    Trace = TraceManager_PSC.Trace(file)


    #Apply a bandpass filter and downsample the data
    T = Trace.Filter(file, pol, drug_arrival, exp_duration, wind_length)
    
    sampling_Hz, X_label, Y_label, time, drug_time, rec_duration, win = T
    
    #Apply a polynomial filter    
    data, global_envlp, precise_envlp, rest = Trace.ModeSubstractor()
    
    #Evaluate the noise unilaterally to establish the detection threshold
    root = Trace.NoiseEvaluator()

    if max(rest) > 3000:
        rest = np.where(rest > 400, 0, rest)
        root = -root
    
    '''                             Analysis                               '''
    height, prominence = 5, root*5
    A = Analysor_PSC_2_0.Analysor(rec_duration, sampling_Hz, 
                             wind_length, drug_arrival)
    
    #Find Events
    threshold, indexes, ampli, amp = A.PeakFinder(rest, root, event_type, pol, height, 
                                      prominence, wind_bl, wind_dr, wind_wsh, wind_length)
    
    #Create a binary list of event and its cumulative function
    peak_binary, spike_tot, ecdf, bl_spikes = A.Binary()
    
    #Calculate the frequency in each bin of time
    (bin_cell_10s, ratio_max_10s, bin_cell_1s)  = A.BinCalculator(bin_10s, bin_1s)

    '''                       Plot section                          '''
    Plot = PlotFactory_PSC_2_0.PlotFactory(rest, data, threshold, indexes,
                                       file_name, exp_duration, rec_duration, 
                                       drug_arrival, wind_length)
    

    #Create an overlay of the curent detected to check their kinetics
    Plot.TracePloter(individual_traces, prominence, height)

    
    '''                           Stock Section                           ''' 
    for k, v in zip(Stock, (file_name, peak_binary, spike_tot, bin_cell_10s, 
                            bin_cell_1s, bl_spikes, ratio_max_10s, 
                            indexes/sampling_Hz, 
                            ecdf, ampli)):
        Stock[k].append(v)
    
    
'''                Create the final Plots                '''
FP = PlotFactory_PSC_2_0.FinalPlots(exp_duration, rec_duration, drug, colorz, 
                                folder_exp_analysis, drug_arrival, event_type, phenotype)

# Show the time course of the frequency variations
FP.TimeCourses(Stock['Bin 10s'], Stock['bl_spikes'], Stock['Ratio 10s'], Stock['Spike Total'], 
               wind_bl/10, wind_dr/10, wind_wsh/10, wind_length/10)

# Create histograms
indiv_val_hz, indiv_val_amp = FP.Histo(Stock['Bin 1s'], Stock['Ampli'],
                     wind_bl, wind_dr, wind_wsh, wind_length, Paired)                                            

df_hz = pd.DataFrame(indiv_val_hz)
writer_exp = pd.ExcelWriter(f'{folder_exp_analysis}/excel_hz.xlsx')
df_hz.to_excel(writer_exp)
writer_exp.save()

df_amp = pd.DataFrame(indiv_val_amp)
writer_exp = pd.ExcelWriter(f'{folder_exp_analysis}/excel_amp.xlsx')
df_amp.to_excel(writer_exp)
writer_exp.save()

winsound.Beep(200, 1000)
