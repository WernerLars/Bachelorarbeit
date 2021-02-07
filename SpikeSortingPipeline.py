import spikeinterface.extractors as se
import spikeinterface.toolkit as st
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw
import json

def createRecordingList():  
    
    return list()

def addRecording(recording_list, recording, sorting_true, study_set_name):  

    recording_name = "recording_"+str(len(recording_list))
    recording_list.append([[recording_name,study_set_name], recording, sorting_true])
    return recording_list

def getSpikeForestStudySetList():
    
    return ["PAIRED_BOYDEN", "PAIRED_CRCNS_HC1", "PAIRED_ENGLISH", "PAIRED_KAMPFF", "PAIRED_MEA64C_YGER", "PAIRED_MONOTRODE",
            "HYBRID_JANELIA", "LONG_DRIFT", "LONG_STATIC", "SYNTH_BIONET", "SYNTH_MAGLAND", "SYNTH_MEAREC_NEURONEXUS",
            "SYNTH_MEAREC_TETRODE", "SYNTH_MONOTRODE", "SYNTH_VISAPY", "MANUAL_FRANKLAB"]

def loadSpikeForestRecordings(recording_list, study_set_name, study_set_number, recording_intervall, download):  
    
    ka.set_config(fr='default_readonly')   
    
    json_file_name = "spikeforest_recordings/recordings/"+study_set_name+"/"+study_set_name+".json"
    json_file = open(json_file_name)
    data = json.load(json_file) 
    
    study = data["studies"][study_set_number]
    recordings = study["recordings"][recording_intervall[0]:recording_intervall[1]]
    
    for rec in recordings:
        recording = AutoRecordingExtractor(rec['directory'], download=download)
        sorting_true = AutoSortingExtractor(rec['firingsTrue'])    
        recording_list = addRecording(recording_list, recording, sorting_true, rec["studyName"]+" : "+rec["name"])   
        
    return recording_list

def printRecordingData(recordings):   
    
    for rec in recordings:  

        print(rec[0][0], ":", rec[0][1])
        recording = rec[1]
        channel_ids = recording.get_channel_ids()
        fs = recording.get_sampling_frequency()
        num_chan = recording.get_num_channels()

        print('Channel ids:', channel_ids)
        print('Sampling frequency:', fs)
        print('Number of channels:', num_chan)
    
        sorting_true = rec[2]
        if sorting_true != "":
            unit_ids = sorting_true.get_unit_ids()
            spike_train = sorting_true.get_unit_spike_train(unit_id=unit_ids[0])

            print('Unit ids:', unit_ids)
            print('Spike train of first unit:', spike_train, "\n")
        
        
# These Methods prints Timeseries,Geometry,Spectrum,Spectogram,
# Raster,Isi Distribution and Correlograms for a list of recordings

def printTimeseries(recordings):
    
    for recording in recordings:  
        plot = sw.plot_timeseries(recording[1])#, channel_ids = [0,1], trange=[0,50])
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"timeseries_"+recording[0][0])
        
def printElectrodeGeometry(recordings):

    for recording in recordings:
        plot = sw.plot_electrode_geometry(recording[1])
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"geometry_"+recording[0][0])
        
def printSpectrum(recordings):
    
    for recording in recordings:
        plot = sw.plot_spectrum(recording[1])
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"spectrum_"+recording[0][0])
        
def printSpectrogram(recordings):
    
    for recording in recordings:
        plot = sw.plot_spectrogram(recording[1], channel=0, nfft=2048)
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"spectrogram_"+recording[0][0])
        
def printRasters(recordings):
    
    for recording in recordings:
        plot = sw.plot_rasters(recording[2], sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"rasters_"+recording[0][0])
        
def printIsiDistribution(recordings):
    
    for recording in recordings:
        plot = sw.plot_isi_distribution(recording[2], bins=10, window=1, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"isi_distribution_"+recording[0][0])
        
def printAutocorrelograms(recordings):

    for recording in recordings:
        plot = sw.plot_autocorrelograms(recording[2], bin_size=1, window=10, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"autocorrelograms_"+recording[0][0])

def printCrosscorrelograms(recordings):
    
    for recording in recordings:
        plot = sw.plot_crosscorrelograms(recording[2], bin_size=0.1, window=5, 
                                  sampling_frequency = recording[1].get_sampling_frequency())
        plot.figure.suptitle(recording[0])
        plot.figure.savefig("figures/"+"crosscorrelograms_"+recording[0][0])
        
  
  
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
            
            if recording[0][0] == key[0]:
                
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
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"unit_waveforms_"+sorter[0]+"_"+rec[0][0][0])
            
            
def printAmplitudeDistribution(sorter_list):

    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            plot = sw.plot_amplitudes_distribution(recording, sorting, max_spikes_per_unit=300)
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"amplitude_distribution_"+sorter[0]+"_"+rec[0][0][0])
    
def printAmplitudeTimeseries(sorter_list):

    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            plot = sw.plot_amplitudes_timeseries(recording, sorting, max_spikes_per_unit=300)
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"amplitude_timeseries_"+sorter[0]+"_"+rec[0][0][0])
    
def printPCAFeatures(sorter_list):
    
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1] 
            plot = sw.plot_pca_features(recording, sorting, colormap='rainbow', nproj=3, max_spikes_per_unit=100)
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"pca_features_"+sorter[0]+"_"+rec[0][0][0])
            
            
def compareWithGroundTruth(sorter_list):
    
    for rec in sorter_list:       
        sorting_true = rec[0][2]
        if not sorting_true == "":
            for sorter in rec[1]:
                sorting = sorter[1]
                comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)
                w_comp = sw.plot_confusion_matrix(comp)
                w_comp.figure.suptitle("GroundTruth_Confusion_Matrix : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
                w_comp.figure.savefig("figures/"+"compare_groundtruth_confusionmatrix_"+sorter[0]+"_"+rec[0][0][0])
                w_agr = sw.plot_agreement_matrix(comp) 
                w_agr.figure.suptitle("GroundTruth_Agreement_Matrix : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
                w_agr.figure.savefig("figures/"+"compare_groundtruth_agreementmatrix_"+sorter[0]+"_"+rec[0][0][0])

def printPerformance(sorter_list):
    
    for rec in sorter_list:
        recording = rec[0][1]
        sorting_true = rec[0][2]
        if not sorting_true == "":
            for sorter in rec[1]:
                sorting = sorter[1]   
                snrs = st.validation.compute_snrs(sorting_true, recording, save_as_property=True)
                comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)               
                w_perf_acc = sw.plot_sorting_performance(comp, property_name='snr', metric='accuracy')
                w_perf_acc.figure.suptitle("Accurcay : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
                w_perf_acc.figure.savefig("figures/"+"performance_accuracy_"+sorter[0]+"_"+rec[0][0][0])               
                w_perf_rec = sw.plot_sorting_performance(comp, property_name='snr', metric='recall')
                w_perf_rec.figure.suptitle("Recall : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
                w_perf_rec.figure.savefig("figures/"+"performance_recall_"+sorter[0]+"_"+rec[0][0][0])              
                w_perf_precision = sw.plot_sorting_performance(comp, property_name='snr', metric='precision')
                w_perf_precision.figure.suptitle("Precision : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0]) 
                w_perf_precision.figure.savefig("figures/"+"performance_precision_"+sorter[0]+"_"+rec[0][0][0])              
                w_perf_false_dis_rate = sw.plot_sorting_performance(comp, property_name='snr', metric='false_discovery_rate')
                w_perf_false_dis_rate.figure.suptitle("False Discovery Rate : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0]) 
                w_perf_false_dis_rate.figure.savefig("figures/"+"performance_false_discovery_rate_"+sorter[0]+"_"+rec[0][0][0])            
                w_perf_miss_rate = sw.plot_sorting_performance(comp, property_name='snr', metric='miss_rate')
                w_perf_miss_rate.figure.suptitle("Miss Rate : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0]) 
                w_perf_miss_rate.figure.savefig("figures/"+"performance_miss_rate_"+sorter[0]+"_"+rec[0][0][0])
            
            
def compareSorters(sorter_entry):
    
    for sorter1 in sorter_entry[1]:
        for sorter2 in sorter_entry[1]:
            if sorter1[0] != sorter2[0]:
                cmp = sc.compare_two_sorters(sorting1=sorter1[1], 
                                             sorting2=sorter2[1],
                                             sorting1_name=sorter1[0],
                                             sorting2_name=sorter2[0])
                plot = sw.plot_agreement_matrix(cmp)
                plot.figure.suptitle(sorter_entry[0][0][0] + " : " + sorter1[0] + " - " + sorter2[0])
                plot.figure.savefig("figures/"+"compare_sorter_"+sorter_entry[0][0][0]+"_"+sorter1[0]+"_"+sorter2[0])
                print(sorter_entry[0][0][0] + " : " + sorter1[0] + " - " + sorter2[0] + "\n")
                print(cmp.match_event_count)
                print(cmp.agreement_scores)
                print("\n")

def compareMultipleSorters(sorter_list):
    
    for rec in sorter_list:
        sorters = list()
        for sorter in rec[1]:
            sorters.append(sorter[1])
        multicomp = sc.compare_multiple_sorters(sorters)
        w_multi = sw.plot_multicomp_graph(multicomp, edge_cmap='coolwarm', node_cmap='viridis', draw_labels=False,
                                  colorbar=True)
        w_multi.figure.suptitle(rec[0][0][0])
        w_multi.figure.savefig("figures/"+"compare_multi_sorters_"+rec[0][0][0])