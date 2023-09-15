# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:19:25 2022

@author: Etienne.CLAUSS
"""


from scipy.signal import find_peaks
import numpy as np


class Analysor():
    def __init__(self, rec_duration, sampling_Hz, wind_length, drug_arrival):
        self.rec_duration = rec_duration
        self.sampling_Hz = sampling_Hz
        self.exp_duration = int(rec_duration/sampling_Hz)
        self.wind_length = wind_length
        self.drug_arrival = drug_arrival    
    
    def PeakFinder(self, rest, root, event_type, pol, prominence, height, 
                    wind_bl, wind_dr, wind_wsh, wind_length, show_fig=False):
        threshold = root*5 if event_type != 'AP' else 10
        self.event_type = event_type
    
        Ind, _ = find_peaks(np.asarray(rest), height = height, prominence=prominence)

        self.indexes = Ind
        
        amp = np.zeros(len(rest))
        amp[:] = np.nan
        for i in self.indexes:
            amp[i] = rest[i]
            
        amp_bl = amp[wind_bl*self.sampling_Hz: wind_bl*self.sampling_Hz + 
                     wind_length*self.sampling_Hz]
        amp_dr = amp[wind_dr*self.sampling_Hz: wind_dr*self.sampling_Hz + 
                     wind_length*self.sampling_Hz]
        amp_wsh = amp[wind_wsh*self.sampling_Hz: wind_wsh*self.sampling_Hz + 
                      wind_length*self.sampling_Hz]
        if not np.nanmean(amp_dr)>0:
            mean_bl = np.nanmean(amp_bl)
            mean_dr = 0
            mean_wsh = np.nanmean(amp_wsh)
            ampli = (mean_bl, mean_dr, mean_wsh)
        else: ampli = (np.nanmean(amp_bl), np.nanmean(amp_dr), np.nanmean(amp_wsh))
        return threshold, self.indexes, ampli, amp 
    

    def Binary(self):
        #Binary list of events
        self.peak_binary = np.zeros(self.rec_duration)
        for index in self.indexes: self.peak_binary[index] += 1
        spike_tot = np.sum(self.peak_binary)
        if spike_tot == 0: spike_tot = 1 
        ecdf = np.cumsum(self.peak_binary)
        ecdf = ecdf/ecdf[-1]
        ecdf = [x for i,x in enumerate(ecdf) if not i%100]
        
        splitted_frame = np.split(self.peak_binary, self.exp_duration)
        mini_bin = [np.nansum(x) for x in splitted_frame]
        bl_spikes = np.sum(mini_bin[:self.drug_arrival]) 

        return self.peak_binary, spike_tot, ecdf, bl_spikes

        
    def BinCalculator(self, bin_10s, bin_1s):
        #Split the data into bins & attribute a binary value to each bin
        #10s bins
        n_bin_10s = int(self.rec_duration/(bin_10s*self.sampling_Hz))
        splitted_cell_binary_10s = np.split(self.peak_binary, n_bin_10s)
        self.bin_cell_10s = [np.nansum(scb) for scb in splitted_cell_binary_10s]  
        ratio_10s = max(self.bin_cell_10s)   
        if ratio_10s == 0: ratio_10s = 1
        #1s bins
        n_bin_1s = int(self.rec_duration/(bin_1s*self.sampling_Hz))
        splitted_cell_binary_1s = np.split(self.peak_binary, n_bin_1s)
        self.bin_cell_1s = [np.nansum(scb) for scb in splitted_cell_binary_1s]  


        return self.bin_cell_10s, ratio_10s, self.bin_cell_1s