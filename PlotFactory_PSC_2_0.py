# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:53:06 2022

@author: Etienne.CLAUSS
"""

import matplotlib.pyplot as plt
import numpy as np, scipy.stats as Stat

class PlotFactory():
    def __init__(self, rest, data, threshold, indexes, file_name, 
                 exp_duration, rec_duration, drug_arrival, wind_length):
        self.rest = rest
        self.data, self.threshold = data, threshold
        self.indexes = indexes
        self.file_name = file_name.split('.')[0]
        self.exp_duration = exp_duration
        self.rec_duration = rec_duration
        self.time = np.linspace(0, self.exp_duration, len(data))
        self.time2 = np.linspace(0, self.exp_duration, self.rec_duration)
        self.drug_arrival = drug_arrival
        self.wind_length = wind_length

    def TracePloter(self, individual_traces, prominence, height, show_fig = False):
        
        box_param = max(self.rest)*0.25
        
        plt.figure().suptitle(self.file_name)
        plt.plot(self.time, self.rest, c='b', lw=.5, label='data')
        plt.axhline(height) 
        plt.text(1000, plt.ylim()[1] + box_param/2, 'Detection parameters')
        plt.text(1000, plt.ylim()[1], f'Prominence: {prominence}')
        plt.text(1000, plt.ylim()[1] - box_param/2, f'Height: {height}')
        plt.axvline(self.drug_arrival, c='gold', lw=2)
        
        for index in self.indexes:
            plt.scatter(self.time[index], self.rest[index], c='r',marker='o')
        plt.savefig(rf'{individual_traces}\{self.file_name}.pdf')
        if not show_fig: plt.close()

        
class FinalPlots():
    def __init__(self, exp_duration, rec_duration, drug, colorz, 
                 folder, drug_arrival, event_type, phenotype):
        self.exp_duration = exp_duration
        self.rec_duration = rec_duration
        self.sampling_Hz = rec_duration/exp_duration
        self.x_ax = np.linspace(0, self.exp_duration, self.rec_duration)
        self.drug = drug
        self.colorz = colorz
        self.folder = folder
        self.drug_arrival = drug_arrival
        self.event_type = event_type
        self.phenotype = phenotype
               
    def TimeCourses(self, bin_10s, bl_spikes, ratio_max, spike_tot, wind_bl, 
                    wind_dr, wind_wsh, wind_length, show_fig = False):
       
        #Time course normalised to the total spikes in the recording
        plt.figure(), plt.title(f'{self.drug}  {self.event_type} {self.phenotype} norm to tot spikes')
        plt.xlabel('Bin 10s'), plt.ylabel('% activity')
        bin_norm = [[(x / spike_tot[y]) for x in bin_10s[y]] 
                    for y, _ in enumerate(bin_10s)]
        bin_norm, x_ax = np.asarray(bin_norm), np.arange(len(bin_norm[0]))
        mean_bin = np.nanmean(bin_norm, axis = 0)
        
        sem_bin = [Stat.sem(bin_norm[:,i]) for i in x_ax]
        plt.fill_between(x_ax, mean_bin-sem_bin, mean_bin+sem_bin)
        plt.axvline(self.drug_arrival/10, c='gold', lw = 2)
        
        plt.axvspan(wind_bl, wind_bl+wind_length,
            color='yellow', alpha=.1)
        plt.axvspan(wind_dr, wind_dr+wind_length,
            color='b', alpha=.1)
        plt.axvspan(wind_wsh, wind_wsh+wind_length,
            color='green', alpha=.1)
        
        plt.plot(x_ax, mean_bin, c='r', zorder=2)      
        plt.savefig(rf'{self.folder}\Time course tot.pdf') 
        if not show_fig: plt.close()
    
        # Individual time course normalised to the total spikes
        for i,x in enumerate(bin_norm):
            plt.figure(), plt.title(f'{self.drug} cell {i+1}')
            plt.xlabel('Bin 10s'), plt.ylabel('% activity')
            plt.axvline(self.drug_arrival/10, c='gold', lw = 2)
            
            plt.axvspan(wind_bl, wind_bl+wind_length,
                color='yellow', alpha=.1)
            plt.axvspan(wind_dr, wind_dr+wind_length,
                color='b', alpha=.1)
            plt.axvspan(wind_wsh, wind_wsh+wind_length,
                color='green', alpha=.1)
            plt.plot(x, c = 'black')
            plt.savefig(rf'{self.folder}\Time course spike tot Freq cell {i+1}.pdf') 
            if not show_fig: plt.close()
           
    
    #Histos
    def Histo(self, bin_1s, ampli, wind_bl, wind_dr, wind_wsh, 
              wind_length, Paired, show_fig = False):        
        x_ax = x_ax = (0, 0.1, 0.2)
        
        #Histo frequency
        bl = [bin_1s[i][wind_bl:wind_bl+wind_length] for i,_ in enumerate(bin_1s)]        
        dr = [bin_1s[i][wind_dr:wind_dr+ wind_length] 
                 for i,_ in enumerate(bin_1s)]        
        wsh = [bin_1s[i][wind_wsh:wind_wsh + wind_length] 
                  for i,_ in enumerate(bin_1s)]        

        bl = [np.mean(bl[i]) for i,_ in enumerate(bl)]
        dr = [np.mean(dr[i]) for i,_ in enumerate(dr)]
        wsh = [np.mean(wsh[i]) for i,_ in enumerate(wsh)]
        val_stat = np.asarray((bl, dr, wsh))  
        indiv_val_hz = val_stat.transpose()
        meanz = np.asarray((np.mean(bl), np.mean(dr), np.mean(wsh)))
        
        bl_sem, dr_sem, wsh_sem = (Stat.sem(bl), Stat.sem(dr), Stat.sem(wsh))  
        semz = np.asarray((np.mean(bl_sem), np.mean(dr_sem), np.mean(wsh_sem)))
        
        plt.figure(), plt.xticks(x_ax, ['Baseline', self.drug, 'Wash']) 
        plt.title(f'{self.drug} {self.event_type} {self.phenotype}')
        plt.ylabel('Frequency (Hz)') 
        plt.bar(x_ax, meanz, width = 0.1,
                yerr=semz, color=self.colorz[1], capsize=10, zorder=0)
        plt.plot(1,0)
        
        a = np.amax(indiv_val_hz)+0.3
        b = (a,a,a)
        plt.plot(x_ax, b, c='k', lw=1, zorder=1)  
        for x in indiv_val_hz:
            [plt.scatter(x_ax[i], x, s=20, c='w', marker='o',
                            edgecolor='k', zorder=2) for i,x in enumerate(x)]  
            plt.plot(x_ax, x, c='k', lw=0.5, zorder=1)   


        plt.savefig(rf'{self.folder}\Histo_wind_fixed.pdf')
        if not show_fig: plt.close()
        
        #Histo amplitude
        bl = [ampli[i][0] for i,_ in enumerate(ampli)]
        dr = [ampli[i][1] for i,_ in enumerate(ampli)]
        wsh = [ampli[i][2] for i,_ in enumerate(ampli)]
        meanz = np.asarray((np.mean(bl), np.mean(dr), np.mean(wsh)))
        
        bl_sem, dr_sem, wsh_sem = (Stat.sem(bl), Stat.sem(dr), Stat.sem(wsh))  
        semz = np.asarray((np.mean(bl_sem), np.mean(dr_sem), np.mean(wsh_sem)))
        
        plt.figure(), plt.xticks(x_ax, ['Baseline', self.drug, 'Wash']) 
        plt.title(f'{self.drug} {self.event_type} {self.phenotype}')
        plt.ylabel('Amplitude (pA)') 
        plt.bar(x_ax, meanz, width = 0.1,
                yerr=semz, color=self.colorz[1], capsize=10, zorder=0)
        plt.plot(1,0)
        
        a = np.amax(ampli)+0.3
        b = (a,a,a)
        plt.plot(x_ax, b, c='k', lw=1, zorder=1)  
        for x in ampli:
            [plt.scatter(x_ax[i], x, s=20, c='w', marker='o',
                            edgecolor='k', zorder=2) for i,x in enumerate(x)]  
            plt.plot(x_ax, x, c='k', lw=0.5, zorder=1)   


        plt.savefig(rf'{self.folder}\Histo_amp_wind_fixed.pdf')
        if not show_fig: plt.close()
        
        return indiv_val_hz, ampli
        