function pli_struct = GENERATE_COHERNECE_MAPS_STRUCT(current_dir,Hd_freq1,Hd_freq2 ,area2load)


data_root           = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
chan_map_folder     = [data_root, '\NPX data\StarndarChannMap'];
area_limit_table    = [data_root, '\Area_limits_GoodLooking.xlsx'];
% npx_raw_data = 
animal_code         = strsplit(current_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};
%% define parameters

sample_length_theta = 4200;
sample_length_delta = 15000;
n_samples           = 500;

%% load synch from synch folder
% load([synch_directory,'\', animal_code, '\synch_model_video2NPX.mat'],'synch_model_video2NPX')




%%  load lfp from current dir
disp('LOADING LFP')
if exist([current_dir,'\','LFP_', area2load, '.mat'], 'file')==2

    NPX_Type        = 2;
    load([current_dir,'\','LFP_', area2load, '.mat'], 'LFP')
elseif exist([current_dir,'\','LFP_', area2load, '.dat'], 'file')==2
    NPX_Type        = 1;
    file_pointer    = fopen([current_dir,'\','LFP_', area2load, '.dat'], 'r');
    LFP             = fread(file_pointer,'int16');
    LFP             = reshape(LFP, 384, numel(LFP)/384);
end


disp('LFP LOADED')
%% select_mid_pag_channel
disp('Loading Channel Map')
areas_by_channel = cell(384,1);
channel_map      = nan(384,2);
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);

% Build animal identifier for area selection
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1
    load([chan_map_folder,'\neuropixPhase3A_kilosortChanMap.mat'], 'xcoords','ycoords', 'chanMap' )
    
 



  for ch_n=1:384
      ch = chanMap(ch_n)+1;
      channel_map(ch,1) = xcoords(ch_n);
      channel_map(ch,2) = ycoords(ch_n);
      areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & ismember(area_limit.Probe_Area, area2load) };
  end

else
    load([current_dir,'\chann_map_PAG.mat'], 'xcoords', 'ycoords','chanMap')

  
    
   for ch_n=1:384
       probe_n = find(any(ismember(hard_coded_x_coords,xcoords(ch_n)),2));
      ch = chanMap(ch_n)+1;
      channel_map(ch,1) = xcoords(ch_n);
      channel_map(ch,2) = ycoords(ch_n);
      areas_by_channel{ch} = area_limit.area{ycoords(ch_n)>=area_limit.depth_start &  ycoords(ch_n)<area_limit.depth_end+1 & area_limit.ProbeNum==probe_n & ismember(area_limit.Probe_Area, area2load)};
  end

      

end
%% obtain_psth



PAG_LFP         = double(LFP);
clear LFP




pli_struct_freq1 = COMPUTE_PLI_WPLI_SEGMENTS(PAG_LFP, Hd_freq1, n_samples, sample_length_delta, [1 5 95 99]);
pli_struct_freq2 = COMPUTE_PLI_WPLI_SEGMENTS(PAG_LFP, Hd_freq2, n_samples, sample_length_theta, [1 5 95 99]);


pli_struct.PLI_matrix_freq1             = pli_struct_freq1.PLI_matrix;  
pli_struct.wPLI_matrix_freq1            = pli_struct_freq1.wPLI_matrix; 
pli_struct.PLI_matrix_pctls_freq1       = pli_struct_freq1.PLI_matrix_pctls; 
pli_struct.wPLI_matrix_pctls_freq1      = pli_struct_freq1.wPLI_matrix_pctls; 
pli_struct.pctls_freq1                  = pli_struct_freq1.pctls;    
pli_struct.freq_range_freq1             = pli_struct_freq1.freq_range;        
pli_struct.Hd_freq1                     = pli_struct_freq1.Hd;   
pli_struct.mean_overlap_freq1           = pli_struct_freq1.mean_overlap;
pli_struct.mean_phaseLag_freq1          = pli_struct_freq1.mean_phaseLag_matrix;
pli_struct.phaseLag_pctls_freq1         = pli_struct_freq1.phaseLag_matrix_pctls;
pli_struct.sample_length_freq1          = sample_length_delta;


pli_struct.PLI_matrix_freq2             = pli_struct_freq2.PLI_matrix;  
pli_struct.wPLI_matrix_freq2            = pli_struct_freq2.wPLI_matrix; 
pli_struct.PLI_matrix_pctls_freq2       = pli_struct_freq2.PLI_matrix_pctls; 
pli_struct.wPLI_matrix_pctls_freq2      = pli_struct_freq2.wPLI_matrix_pctls; 
pli_struct.pctls_freq2                  = pli_struct_freq2.pctls;    
pli_struct.freq_range_freq2             = pli_struct_freq2.freq_range;        
pli_struct.Hd_freq2                     = pli_struct_freq2.Hd; 
pli_struct.mean_overlap_freq2           = pli_struct_freq2.mean_overlap;
pli_struct.mean_phaseLag_freq2          = pli_struct_freq2.mean_phaseLag_matrix;
pli_struct.phaseLag_pctls_freq2         = pli_struct_freq2.phaseLag_matrix_pctls;
pli_struct.sample_length_freq2          = sample_length_theta;




pli_struct.channel_map                  = channel_map;
pli_struct.areas_by_channel             = areas_by_channel;
         

      




end
