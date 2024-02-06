%%
% For Phil's ATG tracking experiemnts
% to tracks previously identified spots


%%

%R_max = 10;

param.mem = 5;
param.dim = 3;
param.good = 5;
param.quiet = 0;


spotData_All = spotData_1_All;



n_frames= length(spotData_All);

spotTableFrames = [];
for ii =1:n_frames
    spotData_1 = spotData_All{ii};
    if ~isempty(spotData_1)
        spotTable = spotData_2_Table(spotData_1);
        spotTableFrames = [spotTableFrames; [ii*ones(size(spotTable(:, 1))), spotTable ]];
    end
end

%XYZT = [spotTableFrames(:, 3:5), spotTableFrames(:, 1)];
%XYZT = spotTableFrames(:, [3:5, 1] );

XYZT = spotTableFrames(:, [3:5, 2, 6:7, 1]);

filepath_XYZT = "\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cropped output\Cell 1 spot data_cropped.xlsx";

writematrix(XYZT, filepath_XYZT);
disp('XYZT file written')

%% if need tracks, run this section

%tracks_data = track(XYZT, R_max, param);

%n_tracks = max(tracks_data(:, end));

%disp(['total tracks: ' num2str(n_tracks) ])

%fieldnames = {'spotData_spotID', 'SI', 'intensityRatio', 'frames'};

%traceList = tracks_2_traceList (tracks_data, 3, fieldnames);


%% Create "Smart average" of LSM data

% INPUT: Provide file and some other input information
%path_1 = 'C:\_Users_\Ivan\Projects\ATG quantification\Test data from Phil\';
path_1 = 'G:\Shared drives\LusKing Shared Drive\Data\Mannino_Phil\LLSM\atg39GFP-nup157Cherry-1\GPUdecon\';

file_type = 'tif'; %'tif'; % or 'dv'
ch_1_token = 'GFP'; % how to recognize filename with oimages form the channel to analyze


% script will read all files in the specified directories, use tokens to choose only desired files
% tokens to exclude filenames
tok_exc_1='.log';
tok_exc_2='TL';
tok_exc_3='._'; % this one to exclude special pseudofiles which refers to moving one level up
tok_exc_4='REF'; %'snapshot'; If don't know what else to exclude put some random set of characters

% tokens to include filenames
tok_inc_1='.tif';
tok_inc_2='ch1';
tok_inc_3='time';


%=========================================================================
% GET FILENAMES in specified folder information

% IMAGE Files selection
filenames = Get_Filenames(path_1);

% use tokens to exclude files
ind_exc_1 = cellfun(@(x) ~isempty(strfind(x,tok_exc_1)),filenames);
ind_exc_2 = cellfun(@(x) ~isempty(strfind(x,tok_exc_2)),filenames);
ind_exc_3 = cellfun(@(x) ~isempty(strfind(x,tok_exc_3)),filenames);
ind_exc_4 = cellfun(@(x) ~isempty(strfind(x,tok_exc_4)),filenames);
%ind6=cellfun(@(x) ~isempty(strfind(x,'Test')),filenames);

% use tokens to include files
ind_inc_1 = cellfun(@(x) ~isempty(strfind(x,tok_inc_1)),filenames);
ind_inc_2 = cellfun(@(x) ~isempty(strfind(x,tok_inc_2)),filenames);
ind_inc_3 = cellfun(@(x) ~isempty(strfind(x,tok_inc_3)),filenames);

% get right names
ind_123 = (~ind_exc_1 & ~ind_exc_2 & ~ind_exc_3 & ~ind_exc_4 & ind_inc_1 & ind_inc_2 & ind_inc_3);
filenames = filenames(ind_123);
filenames_1 = sort(filenames);

%=========================================================================
% LOAD and create 'smart average' in the loop
% to try or debug run by sections whithin the for loop

n_files = length(filenames_1);

Fluo_1_av_TL =[];
for ff=1 :length(filenames_1)
%ff = 1;   

    %========================================================
    %% LOAD IMAGE stack
    
     
    filename_1 = filenames_1{ff};
    %     filename_2=filenames_2{ff};
    %     disp(['Working with ', filename_1])
    
    disp(['Working with ', filename_1])
    
    switch file_type
        case 'tif'
            [Fluo_1, image_info]=load_tiff_stack([path_1,filename_1]);
            %         case 'dv'
            %             % loading using bfread: filename, load 1-st series, all timepoints (or
            %             % specify the range), n-th channel...
            %             channel=ch_1;
            %             Fluo_1 = bfread([path_1,filename_1], 1, 'TimePoints', 'all', 'Channel', channel);
            %             if iscell(Fluo_1)
            %                 Fluo_1=dv_movie_2_matrix(Fluo_1);
            %             end
            %             Fluo_1=uint16(Fluo_1);
    end
        
    disp('DONE image loading')

    im_size = image_info;
    Fluo_1_av = zeros(im_size(1:2));
    for ii = 1: im_size(1)
        for jj = 1:im_size(2)
            Fluo_1_av(ii, jj) = mean(Fluo_1(ii, jj, Fluo_1(ii, jj, :)>0));
            
        end
    end
  
    Fluo_1_av_TL(:, :, ff) = Fluo_1_av; 
    
    
end % enf of files-for-loop

%% SHOW TL-averaged image

im_scale = [0.0, 0.3];
   

%image  = mean(Fluo_1_av_TL, 3);

image_1  = Fluo_1_av_TL (:, :, 45);


im_min = min(image_1(:));
    im_max = max(image_1(:));
    im_range = im_max-im_min;
    im_0 = im_min+im_scale(1)*im_range;
    im_1 = im_min+im_scale(2)*im_range;
   
    figure;
        imshow(image_1,[im_0, im_1],'InitialMagnification',200,'Border','tight');
        hold on
  
    
    


