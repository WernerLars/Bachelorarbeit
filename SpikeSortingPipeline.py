from spikeforest2_utils import AutoRecordingExtractor, AutoSortingExtractor
import kachery as ka
import spikeinterface.extractors as se
import spikeinterface.toolkit as st
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw
import json

def createRecordingList():  
    
    # No Parameter 
    # Returns a list for storing Recordings. You can also just define a list directly without using this method.
    
    return list()

def addRecording(recording_list, recording, sorting_true, study_set_name):  

    # Parameter:
        # recording_list    : this is the data object, where the recording is stored
        # recording         : the RecordingExtractor, that we want to store 
        # sorting_true      : the SortingExtractor, that if available we want to store
        # study_set_name    : for storing the study set name
    # Returns the recording_list with the stored elements.

    # elements in the list are counted and integrated in the recording_name, will be used to compare names in createSorterList
    recording_name = "recording_"+str(len(recording_list))   

    # Appends the elements to the list as an array (recording_name and study_set_name as an extra array)
    recording_list.append([[recording_name,study_set_name], recording, sorting_true])
       
    return recording_list

def getSpikeForestStudySetList():
    
    # No Parameter
    # Returns all available study_set_names from spikeforest2 as a list
    
    return ["PAIRED_BOYDEN", "PAIRED_CRCNS_HC1", "PAIRED_ENGLISH", "PAIRED_KAMPFF", "PAIRED_MEA64C_YGER", "PAIRED_MONOTRODE",
            "HYBRID_JANELIA", "LONG_DRIFT", "LONG_STATIC", "SYNTH_BIONET", "SYNTH_MAGLAND", "SYNTH_MEAREC_NEURONEXUS",
            "SYNTH_MEAREC_TETRODE", "SYNTH_MONOTRODE", "SYNTH_VISAPY", "MANUAL_FRANKLAB"]

def loadSpikeForestRecordings(recording_list, study_set_name, study_set_number, recording_intervall, download):  
    
    # Parameter:
        # recording_list        : this is the data object, where the recordings are stored
        # study_set_name        : to load the right path and json file
        # study_set_number      : to get the study in the study_set you want the recordings from
        # recording_intervall   : to have a boundary for slicing recordings
        # download              : boolean if you want to store the recording files on your computer
    # Returns the recording_list with all new added recordings from spikeforest2
    
    # From SpikeForest2
    ka.set_config(fr='default_readonly')   
    
    # Loading the json File (dictionary) from spikeforest_recordings
    json_file_name = "spikeforest_recordings/recordings/"+study_set_name+"/"+study_set_name+".json"
    json_file = open(json_file_name)
    data = json.load(json_file) 
    
    # Extracting the recordings from the right study (look up the structure of the json files for understanding this process)
    study = data["studies"][study_set_number]
    recordings = study["recordings"][recording_intervall[0]:recording_intervall[1]]  # Slicing the study
    
    # Every Recording as SHA-1 URIs for the recording and the ground truth data. Loading the data with the Extractors from spikeforest2
    # and adding every recording to the recording_list. Look in JSON File for the structure of the recording dictionary.
    for rec in recordings:
        recording = AutoRecordingExtractor(rec['directory'], download=download) 
        sorting_true = AutoSortingExtractor(rec['firingsTrue'])    
        recording_list = addRecording(recording_list, recording, sorting_true, rec["studyName"]+" : "+rec["name"])   
        
    return recording_list

def printRecordingData(recordings):   
    
    # Parameter
        # recordings  : a list with stored recordings
    # No Returns

    # For every entry in the list the name, channel ids, frequency and number of channels from the recording are printed  
    for rec in recordings:  

        print(rec[0][0], ":", rec[0][1])
        recording = rec[1]
        channel_ids = recording.get_channel_ids()
        fs = recording.get_sampling_frequency()
        num_chan = recording.get_num_channels()

        print('Channel ids:', channel_ids)
        print('Sampling frequency:', fs)
        print('Number of channels:', num_chan)
 
     # If we have ground truth, then unit ids and the first spike train will be printed
        sorting_true = rec[2]
        if sorting_true != "":
            unit_ids = sorting_true.get_unit_ids()
            spike_train = sorting_true.get_unit_spike_train(unit_id=unit_ids[0])

            print('Unit ids:', unit_ids)
            print('Spike train of first unit:', spike_train, "\n")

def printTimeseries(recordings):
    
    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Timeseries for every recording in the list.
   
    for recording in recordings:  
        plot = sw.plot_timeseries(recording[1],  channel_ids = [0,1] , trange=[0,0.01])
        
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"timeseries_"+recording[0][0])
  
def printElectrodeGeometry(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Electrode Geometry for every recording in the list.
   
    for recording in recordings:
        plot = sw.plot_electrode_geometry(recording[1])
        
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"geometry_"+recording[0][0])

def printSpectrum(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Spectrum for every recording in the list.
   
    for recording in recordings:
        plot = sw.plot_spectrum(recording[1])
        
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"spectrum_"+recording[0][0])

def printSpectrogram(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Spectrogram for every recording in the list.
   
    
    for recording in recordings:
        plot = sw.plot_spectrogram(recording[1], channel=0, nfft=2048)
        
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"spectrogram_"+recording[0][0])

def printRasters(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Rasters for every recording in the list.
   
    
    for recording in recordings:
        plot = sw.plot_rasters(recording[2], sampling_frequency = recording[1].get_sampling_frequency())
        
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"rasters_"+recording[0][0])

def printIsiDistribution(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the IsiDistribution for every recording in the list.
   
    
    for recording in recordings:
        plot = sw.plot_isi_distribution(recording[2], bins=10, window=1, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
                                 
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"isi_distribution_"+recording[0][0])

def printAutocorrelograms(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Autocorrelograms for every recording in the list.
   

    for recording in recordings:
        plot = sw.plot_autocorrelograms(recording[2], bin_size=1, window=10, 
                                 sampling_frequency = recording[1].get_sampling_frequency())
                                 
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0][0]+" : "+recording[0][1])
        plot.figure.savefig("figures/"+"autocorrelograms_"+recording[0][0])

def printCrosscorrelograms(recordings):

    # Parameter
        # recordings : a list with stored recordings
    # No Returns
    
    # Method for creating the graph of the Crosscorrelograms for every recording in the list.
   
    
    for recording in recordings:
        plot = sw.plot_crosscorrelograms(recording[2], bin_size=0.1, window=5, 
                                  sampling_frequency = recording[1].get_sampling_frequency())
                                  
        # Using the name of the recording and the study for naming the pictures and files.
        plot.figure.suptitle(recording[0])
        plot.figure.savefig("figures/"+"crosscorrelograms_"+recording[0][0])

def runSpikeSorting(recordings,working_folder):
    
    # Parameter
        # recordings     : list of recordings
        # working_folder : string for the name of the folder, where the sorting results are stored (folder must not exist)
    # Returns the spike sorting results as a dictionary
    
    
    # In recordings, we have also information about the name and SortingExtractors of the recordings. 
    # The run_sorters method from spikeinterface only allow RecordingExtractors, so we must filter out the RecordingExtractors of the list.
    recording_list = list()
    for recording in recordings:
        recording_list.append(recording[1])
    
    
    # With this code below, spike sorting will be executed in spikeinterface with all installed sorters and the filtered Recordinglist
    spike_sorting = ss.run_sorters(sorter_list = ss.installed_sorters(),
                               recording_dict_or_list = recording_list,
                               working_folder = working_folder)
    
    return spike_sorting

def createSorterList(recordings,spike_sorting):
    
    # Parameter
        # recordings        : a list of recordings
        # spike_sorting     : the results of the runSpikeSorting method
    # Returns a list, where recordings and spike_sorting will be combined, so all data is in one list
    
    
    # New list, in which the other two lists will be combined
    sorter_list = list()
    
    # Every recording in recordings will be combined with all matching results of the spike_sorting list
    for recording in recordings:
        
        # New List where the SortingExtractors of the recording and their name will be stored.
        sorters = list()
        
        # We have to check all SortingExtractors, if it is a result for the recording
        for key in spike_sorting:
            
            # To match a recording with its result, there has to be a key in both lists and therefore 
            # the name of the recording was choosen, because it is in both lists (thats why the name was computed in addRecording).
            if recording[0][0] == key[0]:
                
                # Storing name and SortingExtractor as an array
                sorters.append([key[1],spike_sorting[key]])
        
        # After all SortingExtractors of an recording has been found, it will be stored into the sorter_list
        sorter_list.append([recording,sorters])
        
    return sorter_list

def printSorterList(sorter_list):
    
    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    
    # Prints every entry of the list
    for sorter in sorter_list:
        print(sorter,"\n")

def printSorterData(sorter_list):
           
    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # For every recording of the list we compute and print the 
    # signal-to-noise ratios, isi_violations and isolation_distances for every sorter result to its recording
    for entry in sorter_list:

        print("-"*50,"\n")
        print(entry[0][0],"\n")
        print("-"*50,"\n")
        
        for sorter in entry[1]:
            
            print(sorter[0],"\n")
            sorting = sorter[1]
            recording = entry[0][1]
            
            # Here the metrics are computed for the specific spike sorting result
            snrs = st.validation.compute_snrs(sorting, recording)
            isi_violations = st.validation.compute_isi_violations(sorting, 
                                                          duration_in_frames=recording.get_num_channels())
            isolations = st.validation.compute_isolation_distances(sorting, recording)
            
            # Here the metrics are printed
            print('SNR', snrs,"\n")
            print('ISI violation ratios', isi_violations,"\n")
            print('Isolation distances', isolations, "\n"*2)

def printUnitWaveforms(sorter_list):
    
    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # Prints for every sorting result the unit waveforms
    for rec in sorter_list:       
        recording = rec[0][1]    
        
        for sorter in rec[1]:
            sorting = sorter[1] 
            plot = sw.plot_unit_waveforms(recording, sorting, max_spikes_per_unit=100)
            
            # Using the name of the recording and the sorter for naming pictures and files.
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"unit_waveforms_"+sorter[0]+"_"+rec[0][0][0])
  
def printAmplitudeDistribution(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # Prints for every sorting result the amplitude distribution
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            plot = sw.plot_amplitudes_distribution(recording, sorting, max_spikes_per_unit=300)
            
            # Using the name of the recording and the sorter for naming pictures and files.           
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"amplitude_distribution_"+sorter[0]+"_"+rec[0][0][0])
    
def printAmplitudeTimeseries(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns

    # Prints for every sorting result the amplitude timeseries
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1]   
            
            plot = sw.plot_amplitudes_timeseries(recording, sorting, max_spikes_per_unit=300)
            
            # Using the name of the recording and the sorter for naming pictures and files.
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"amplitude_timeseries_"+sorter[0]+"_"+rec[0][0][0])
    
def printPCAFeatures(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # Prints for every sorting result the pca features
    for rec in sorter_list:       
        recording = rec[0][1]        
        for sorter in rec[1]:
            sorting = sorter[1] 
            
            plot = sw.plot_pca_features(recording, sorting, colormap='rainbow', nproj=3, max_spikes_per_unit=100)
            
            # Using the name of the recording and the sorter for naming pictures and files.
            plot.figure.suptitle(rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])
            plot.figure.savefig("figures/"+"pca_features_"+sorter[0]+"_"+rec[0][0][0])
     
def compareWithGroundTruth(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    for rec in sorter_list:       
        sorting_true = rec[0][2]
        
        # If we dont have ground truth, the sorting result of the sorter will be ignored
        if not sorting_true == "":
            for sorter in rec[1]:
                sorting = sorter[1]
                
                # for comparing the sorting results with ground truth
                comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)
                
                # Plotting Confusion Matrix
                w_comp = sw.plot_confusion_matrix(comp)
                
                # Using the name of the recording and the sorter for naming pictures and files.
                w_comp.figure.suptitle("GroundTruth_Confusion_Matrix : "+rec[0][0][0]+" : "+sorter[0])
                w_comp.figure.savefig("figures/"+"compare_groundtruth_confusionmatrix_"+sorter[0]+"_"+rec[0][0][0])
                
                # Plotting Agreement Matrix
                w_agr = sw.plot_agreement_matrix(comp) 
                
                # Using the name of the recording and the sorter for naming pictures and files.
                w_agr.figure.suptitle("GroundTruth_Agreement_Matrix : "+rec[0][0][0]+" : "+sorter[0])
                w_agr.figure.savefig("figures/"+"compare_groundtruth_agreementmatrix_"+sorter[0]+"_"+rec[0][0][0])

def printPerformance(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # Printing metrics with ground truth comparison for every sorter that has ground truth data
    for rec in sorter_list:
        recording = rec[0][1]
        sorting_true = rec[0][2]
        
        # If we dont have ground truth, the sorting result of the sorter will be ignored
        if not sorting_true == "":
            for sorter in rec[1]:
                sorting = sorter[1]   
                
                # Computing snrs as property and comp for comparing the sorting results with ground truth
                snrs = st.validation.compute_snrs(sorting_true, recording, save_as_property=True)
                comp = sc.compare_sorter_to_ground_truth(sorting_true,sorting)       
                
                # Metric = Accurcay
                w_perf_acc = sw.plot_sorting_performance(comp, property_name='snr', metric='accuracy')
                # Using the name of the recording and the sorter for naming pictures and files.
                w_perf_acc.figure.suptitle("Accurcay : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])                
                w_perf_acc.figure.savefig("figures/"+"performance_accuracy_"+sorter[0]+"_"+rec[0][0][0])  
                
                # Metric = Recall
                w_perf_rec = sw.plot_sorting_performance(comp, property_name='snr', metric='recall')
                # Using the name of the recording and the sorter for naming pictures and files.
                w_perf_rec.figure.suptitle("Recall : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])                
                w_perf_rec.figure.savefig("figures/"+"performance_recall_"+sorter[0]+"_"+rec[0][0][0]) 
                
                # Metric = Precision
                w_perf_precision = sw.plot_sorting_performance(comp, property_name='snr', metric='precision')
                # Using the name of the recording and the sorter for naming pictures and files.
                w_perf_precision.figure.suptitle("Precision : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0]) 
                w_perf_precision.figure.savefig("figures/"+"performance_precision_"+sorter[0]+"_"+rec[0][0][0])   

                # Metric = False Discovery Rate
                w_perf_false_dis_rate = sw.plot_sorting_performance(comp, property_name='snr', metric='false_discovery_rate')
                # Using the name of the recording and the sorter for naming pictures and files.
                w_perf_false_dis_rate.figure.suptitle("False Discovery Rate : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0]) 
                w_perf_false_dis_rate.figure.savefig("figures/"+"performance_false_discovery_rate_"+sorter[0]+"_"+rec[0][0][0])  

                # Metric = Miss Rate
                w_perf_miss_rate = sw.plot_sorting_performance(comp, property_name='snr', metric='miss_rate')
                # Using the name of the recording and the sorter for naming pictures and files.
                w_perf_miss_rate.figure.suptitle("Miss Rate : "+rec[0][0][0]+" : "+rec[0][0][1]+" : "+sorter[0])                 
                w_perf_miss_rate.figure.savefig("figures/"+"performance_miss_rate_"+sorter[0]+"_"+rec[0][0][0])
                        
def compareSorters(sorter_entry):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns
    
    # Comparing every sorter with each other
    for sorter1 in sorter_entry[1]:
        for sorter2 in sorter_entry[1]:
        
            # To ignore comparison between the same sorter
            if sorter1[0] != sorter2[0]:
                
                # Two sorters will be compared in spikeinterface
                cmp = sc.compare_two_sorters(sorting1=sorter1[1], 
                                             sorting2=sorter2[1],
                                             sorting1_name=sorter1[0],
                                             sorting2_name=sorter2[0])
                                             
                # Creating Agreement Matrix
                plot = sw.plot_agreement_matrix(cmp)
                plot.figure.suptitle(sorter_entry[0][0][0] + " : " + sorter1[0] + " - " + sorter2[0])
                plot.figure.savefig("figures/"+"compare_sorter_"+sorter_entry[0][0][0]+"_"+sorter1[0]+"_"+sorter2[0])
                
                # Print match event counts and agreement_scores of the two sorters
                print(sorter_entry[0][0][0] + " : " + sorter1[0] + " - " + sorter2[0] + "\n")
                print(cmp.match_event_count)
                print(cmp.agreement_scores)
                print("\n")

def compareMultipleSorters(sorter_list):

    # Parameter
        # sorter_list   : list of recording with its results
    # No Returns

    # Printing Multicomparion graphs
    for rec in sorter_list:
        sorters = list()
        for sorter in rec[1]:
            sorters.append(sorter[1])
        multicomp = sc.compare_multiple_sorters(sorters)
        w_multi = sw.plot_multicomp_graph(multicomp, edge_cmap='coolwarm', node_cmap='viridis', draw_labels=False,
                                  colorbar=True)
        w_multi.figure.suptitle(rec[0][0][0])
        w_multi.figure.savefig("figures/"+"compare_multi_sorters_"+rec[0][0][0])