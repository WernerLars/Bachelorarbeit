import spikeinterface.extractors as se
import spikeinterface.toolkit as st
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw
import json

def createRecordingList():  
    
    # Returns a List.
    
    return list()

def addRecording(recordings,recording_name,recording,sorting_true):
    
    # Input: 
        #recordings is a list
        #recording_name is a string, must be in Format "recording_x" with x a number. 
        #                If you add recordings count up from 0
        #recording is an RecordingExtractor
        #sorting_true is Ground Truth Data
        
    # Output is the recording list with added recording
    
    recordings.append([recording_name,recording,sorting_true])
    return recordings

def loadSpikeForestRecordings(recordings,json_file_path,json_file_name,download,number_of_rec):  
    
    # Input
        #recordings 
        #json_file_path in the spikeforest_recordings folder
        #json_file_name is the name of the json file
        #download  True or False, depends if you want to download the recordings
        #number_of_rec is the number of recordings you want to get of the json file
    
    # loads the SpikeForestRecordings registered in the json file from top to bottom
    
    # Output is the recording list with all SpikeForest Recordings you selected
    
    ka.set_config(fr='default_readonly')   
    json_file = open(json_file_path+json_file_name)
    data = json.load(json_file)  
    
    recording_name = "recording_"
    count_recording = len(recordings)
    for rec in data["recordings"][:number_of_rec]:
        recording = AutoRecordingExtractor(rec['directory'],download=download)
        sorting_true = AutoSortingExtractor(rec['firingsTrue'])
        recordings = addRecording(recordings,recording_name+str(count_recording),recording,sorting_true)   
        count_recording += 1
    return recordings

def printRecordingData(recordings):   
    
    # Input is an recording list
    
    # Prints Channel Id, Sampling Frequency and the Number of Channel of all Recordings registered
    
    for rec in recordings:  
        
        print(rec[0])
        recording = rec[1]
        channel_ids = recording.get_channel_ids()
        fs = recording.get_sampling_frequency()
        num_chan = recording.get_num_channels()

        print('Channel ids:', channel_ids)
        print('Sampling frequency:', fs)
        print('Number of channels:', num_chan)
    
        sorting_true = rec[2]
        unit_ids = sorting_true.get_unit_ids()
        spike_train = sorting_true.get_unit_spike_train(unit_id=unit_ids[0])

        print('Unit ids:', unit_ids)
        print('Spike train of first unit:', spike_train, "\n")
        
        
# These Methods prints Timeseries,Geometry,Spectrum,Spectogram,
# Raster,Isi Distribution and Correlograms for a list of recordings



def printTimeseries(recordings):
    
    for recording in recordings:  
        plot = sw.plot_timeseries(recording[1], trange=[0,5])
        plot.figure.suptitle(recording[0])
        
def printElectrodeGeometry(recordings):

    for recording in recordings:
        plot = sw.plot_electrode_geometry(recording[1])
        plot.figure.suptitle(recording[0])
        
def printSpectrum(recordings):
    
    for recording in recordings:
        plot = sw.plot_spectrum(recording[1])
        plot.figure.suptitle(recording[0])
        
def printSpectrogram(recordings):
    
    for recording in recordings:
        plot = sw.plot_spectrogram(recording[1], channel=0, nfft=2048)
        plot.figure.suptitle(recording[0])
        
def printRasters(recordings):
    
    for recording in recordings:
        plot = sw.plot_rasters(recording[2], sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0])
        
def printIsiDistribution(recordings):
    
    for recording in recordings:
        plot = sw.plot_isi_distribution(recording[2], bins=10, window=1, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0])
        
def printAutocorrelograms(recordings):

    for recording in recordings:
        plot = sw.plot_autocorrelograms(recording[2], bin_size=1, window=10, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0])

def printCrosscorrelograms(recordings):
    
    for recording in recordings:
        plot = sw.plot_crosscorrelograms(recording[2], bin_size=0.1, window=5, 
                                  sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0])
        
  

  
# These methods are for running the Spike Sorting for a list of recordings with all installed Spike Sorters,
# creating a list that contains all recordings and their sorters with sorting results
# You can also print the SorterList and some SorterData



def runSpikeSorting(recordings,working_folder):
    
    recording_list = list()
    for recording in recordings:
        recording_list.append(recording[1])
    
    spike_sorting = ss.run_sorters(sorter_list = ss.installed_sorters(),
                               recording_dict_or_list = recording_list,
                               working_folder = working_folder)
    
    return spike_sorting

def createSorterList(recordings,spike_sorting):
    
    sorter_list = list()
    
    for recording in recordings:
        
        sorters = list()
        
        for key in spike_sorting:
            
            if recording[0] == key[0]:
                
                sorters.append([key[1],spike_sorting[key]])
            
        sorter_list.append([recording,sorters])
        
    return sorter_list

def printSorterList(sorter_list):
    
    for sorter in sorter_list:
        print(sorter,"\n")
        
def printSorterData(sorter_list):
    
    for entry in sorter_list:

        print("-"*50,"\n")
        print(entry[0][0],"\n")
        print("-"*50,"\n")
        
        for sorter in entry[1]:
            
            print(sorter[0],"\n")
            sorting = sorter[1]
            recording = entry[0][1]
    
            snrs = st.validation.compute_snrs(sorting, recording)
            isi_violations = st.validation.compute_isi_violations(sorting, 
                                                          duration_in_frames=recording.get_num_channels())
            isolations = st.validation.compute_isolation_distances(sorting, recording)

            print('SNR', snrs,"\n")
            print('ISI violation ratios', isi_violations,"\n")
            print('Isolation distances', isolations, "\n"*2)
     

     
# These are Functions to print Graphs of the Unit Waveforms, Amplitude Distribution,
# Amplitude Timeseries and PCA Features for a sorterlist



def printUnitWaveforms(sorter_list):
    
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1] 
            plot = sw.plot_unit_waveforms(recording, sorting, max_spikes_per_unit=100)
            plot.figure.suptitle(rec[0][0]+" : "+sorter[0])
            
def printAmplitudeDistribution(sorter_list):

    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            plot = sw.plot_amplitudes_distribution(recording, sorting, max_spikes_per_unit=300)
            plot.figure.suptitle(rec[0][0]+" : "+sorter[0])
    
def printAmplitudeTimeseries(sorter_list):

    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            plot = sw.plot_amplitudes_timeseries(recording, sorting, max_spikes_per_unit=300)
            plot.figure.suptitle(rec[0][0]+" : "+sorter[0])
    
def printPCAFeatures(sorter_list):
    
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1] 
            plot = sw.plot_pca_features(recording, sorting, colormap='rainbow', nproj=3, max_spikes_per_unit=100)
            plot.figure.suptitle(rec[0][0]+" : "+sorter[0])
            
            
def compareWithGroundTruth(sorter_list):
    
    for rec in sorter_list:       
        sorting_true = rec[0][2]
        for sorter in rec[1]:
            sorting = sorter[1]
            comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)
            w_comp = sw.plot_confusion_matrix(comp)
            w_comp.figure.suptitle(rec[0][0] + " : " + sorter[0] + " - Confusion Matrix")
            w_agr = sw.plot_agreement_matrix(comp) 
            w_agr.figure.suptitle(rec[0][0] + " : " + sorter[0] + " - Agreement Matrix")

def printPerformance(sorter_list):
    
    for rec in sorter_list:
        recording = rec[0][1]
        sorting_true = rec[0][2]
        for sorter in rec[1]:
            sorting = sorter[1]   
            snrs = st.validation.compute_snrs(sorting_true, recording, save_as_property=True)
            comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)
            w_perf_acc = sw.plot_sorting_performance(comp, property_name='snr', metric='accuracy')
            w_perf_acc.figure.suptitle(rec[0][0] + " : " + sorter[0])
            w_perf_rec = sw.plot_sorting_performance(comp, property_name='snr', metric='recall')
            w_perf_rec.figure.suptitle(rec[0][0] + " : " + sorter[0])
            w_perf_precision = sw.plot_sorting_performance(comp, property_name='snr', metric='precision')
            w_perf_precision.figure.suptitle(rec[0][0] + " : " + sorter[0])  
            
def compareSorters(rec):
    
    for sorter1 in rec[1]:
        for sorter2 in rec[1]:
            if sorter1[0] != sorter2[0]:
                cmp = sc.compare_two_sorters(sorting1=sorter1[1], 
                                             sorting2=sorter2[1],
                                             sorting1_name=sorter1[0],
                                             sorting2_name=sorter2[0])
                plot = sw.plot_agreement_matrix(cmp)
                plot.figure.suptitle(rec[0][0] + " : " + sorter1[0] + " - " + sorter2[0])
                print(cmp.match_event_count)
                print(cmp.agreement_scores)

def compareMultipleSorters(sorter_list):
    
    for rec in sorter_list:
        sorters = list()
        for sorter in rec[1]:
            sorters.append(sorter[1])
        multicomp = sc.compare_multiple_sorters(sorters)
        w_multi = sw.plot_multicomp_graph(multicomp, edge_cmap='coolwarm', node_cmap='viridis', draw_labels=False,
                                  colorbar=True)
        w_multi.figure.suptitle(rec[0][0])