
% Script to detect spots regardless of their shape using 'findIrregularSpots3'
% goes file-by-file in the loop,; wihtin the loop:
%  - load image file
% - detect spots according to parameters
% - merge proximal spots (if asked)
% - save spots data into hierarchical spotData_all array


%% INTIALIZATION: specify general parameters, constants etc
 
% pixel and dz size, not reuired for the spot detection but will be saved with the results
dx=0.1; % pixel size
dz=0.2; % z-stack step

%=========================================================================

%% INPUT: Provide file and some other input information

Fluo_offset = 10;

path_1 = 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\Colocalization stills\fis1KO DNM1WT- 04042023\Cell 2\Redo spot finder\mCherry\Input\';

file_type = 'tif'; %'tif'; % or 'dv'
%ch_1_token = 'ch1'; % how to recognize filename with oimages form the channel to analyze

% filtering, merging etc
filter_spots_1 = 0; % whether to filter Spots by Intensity and IntensityRatio. If so, provide low_cutoffs
%filter_spots_2=0; % whether to filter Spots by Intensity and IntensityRatio. If so, provide low_cutoffs
% cutoffs:
SI_cutoff = 400;
IR_cutoff = 1.7;
spot_merge = 0; % whthere to merge spots that are too close to each other
R_merge = 5; % if spots a re closer than R_max they will be merged

% location where to save otuputs
ssave = 1;
path_2 = 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\Colocalization stills\fis1KO DNM1WT- 04042023\Cell 2\Redo spot finder\mCherry\Output\'; %path to save spot detection results 
path_3  = 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\Colocalization stills\fis1KO DNM1WT- 04042023\Cell 2\Redo spot finder\mCherry\Output\'; % path to save max-projectin images with spots detected
%path_3=[path_1,'/Spots'];
filename_out='Cell 7 tp 41 and 42 mCh foci.mat';

% whether to show or not images etc
show_fig = 1; % to show or not rsults of spot detection for each
show_pars=1; % to show or not spotIntensity (top-right)  and intensityRatio (top-left) from spotList
% for spots visualization
im_sc=[0, 0.5]; % image scaling as [black_px_intenisty, white_px_intensity]

save_fig = 1;
tag = ''; % to tag some saved figures

% script will read all files in the specified directories, use tokens to choose only desired files
% tokens to exclude filenames
tok_exc_1='.log';
tok_exc_2='TL';
tok_exc_3='._'; % this one to exclude special pseudofiles which refers to moving one level up
tok_exc_4='REF'; %'snapshot'; If don't know what else to exclude put some random set of characters

% tokens to include filenames
tok_inc_1='.tif';
tok_inc_2='ch1';
tok_inc_3='488nm';


%=========================================================================

%% GET FILENAMES in specified folder information

% IMAGE Files selection
filenames = Get_Filenames(path_1);

% use tokens to exclude files
%ind_exc_1 = cellfun(@(x) ~isempty(strfind(x,tok_exc_1)),filenames);
%ind_exc_2 = cellfun(@(x) ~isempty(strfind(x,tok_exc_2)),filenames);
%ind_exc_3 = cellfun(@(x) ~isempty(strfind(x,tok_exc_3)),filenames);
%ind_exc_4 = cellfun(@(x) ~isempty(strfind(x,tok_exc_4)),filenames);
%ind6=cellfun(@(x) ~isempty(strfind(x,'Test')),filenames);

% use tokens to include files
%ind_inc_1 = cellfun(@(x) ~isempty(strfind(x,tok_inc_1)),filenames);
%ind_inc_2 = cellfun(@(x) ~isempty(strfind(x,tok_inc_2)),filenames);
%ind_inc_3 = cellfun(@(x) ~isempty(strfind(x,tok_inc_3)),filenames);

% get right names
%ind_123 = (~ind_exc_1 & ~ind_exc_2 & ~ind_exc_3 & ~ind_exc_4 & ind_inc_1 & ind_inc_2 & ind_inc_3);
%filenames = filenames(ind_123);

filenames_1 = sort(filenames);

%=========================================================================
%% LOAD and detect SPOTS in the loop
% to try or debug run by sections whithin the for loop

spotData_1_All={};
%spotData_2_All={};
removed_spots=[];

for ff=1 :length(filenames_1)
    
    %========================================================
    %% LOAD IMAGE stack
    filename_1=filenames_1{ff};
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
    
    Fluo_1 = Fluo_1 + Fluo_offset;


    
%     switch file_type
%         case 'tif'
%             [Fluo_2, image_info]=load_tiff_stack([path_1,filename_2]);
%             %         case 'dv'
%             %             % loading using bfread: filename, load 1-st series, all timepoints (or
%             %             % specify the range), n-th channel...
%             %             channel=ch_2;
%             %             Fluo_2 = bfread([path_1,'/',filename_2], 1, 'TimePoints', 'all', 'Channel', channel);
%             %             if iscell(Fluo_2)
%             %                 Fluo_2=dv_movie_2_matrix(Fluo_2);
%             %             end
%             %             Fluo_2=uint16(Fluo_2);
%     end
%     
    disp('DONE image loading')
    
    
    %========================================================
    %% LOAD or DETECT SPOTS
    %Provide  para meters for spot detection
    param.peakRadius=3.5;
    param.intensityRatioThreshold=1.5;
    param.shellThickness=1;
    param.edgeDist=7;
    param.centerDist=3;
    param.fitRadius=1.0;

    
    % spot detection itself
    disp('Starting spot detection...')
    tic
    [spotData_1,spotDetection_params_1]= findIrregularSpots3(Fluo_1,'peakRadius',param.peakRadius,...
        'intensityRatioThreshold',param.intensityRatioThreshold,...
        'shellThickness',param.shellThickness,...
        'edgeDist',param.edgeDist,...
        'centerDist',param.centerDist, ...
        'fitRadius',param.fitRadius);
    
    %         [spotData_2,spotDetection_params_2]= findIrregularSpots3(Fluo_2,'peakRadius',param.peakRadius,...
    %             'intensityRatioThreshold',param.intensityRatioThreshold,...
    %             'shellThickness',param.shellThickness,...
    %             'edgeDist',param.edgeDist,...
    %             'centerDist',param.centerDist, ...
    %             'fitRadius',param.fitRadius);
    
    % collect results
    spotData_1_All{ff}=spotData_1;
    %         spotData_2_All{ff}=spotData_2;
    disp('... DONE with the spot detection!!!')
    toc
    
    beep
    
    %========================================================
    %% FILTER SPOTS (if asked) by Intensity or Ratio
    if filter_spots_1
        spotData_F={};
        ss=0;
        for jj=1:length(spotData_1)
            if spotData_1{jj}.spotIntensity>SI_cutoff && spotData_1{jj}.intensityRatio>IR_cutoff
                ss=ss+1;
                spotData_F{ss}=spotData_1{jj};
            end
            
        end
        % keep track of filter pars
        filter_pars.SI_cutoff=SI_cutoff;
        filter_pars.IR_cutoff=IR_cutoff;
        % keep old spotData, JIC
        spotData_0=spotData_1;
        spotData_1=spotData_F;
    else
        filter_pars = [];
    end
    filter_pars_1 = filter_pars;
    
    %     if filter_spots_2
    %         spotData_F={};
    %         ss=0;
    %         for jj=1:length(spotData_2)
    %             if spotData_2{jj}.spotIntensity>SI_cutoff && spotData_2{jj}.intensityRatio>IR_cutoff
    %                 ss=ss+1;
    %                 spotData_F{ss}=spotData_2{jj};
    %             end
    %
    %         end
    %         % keep track of filter pars
    %         filter_pars.SI_cutoff=SI_cutoff;
    %         filter_pars.IR_cutoff=IR_cutoff;
    %         % keep old spotData, JIC
    %         spotData_0=spotData_2;
    %         spotData_2=spotData_F;
    %     else
    %         filter_pars = [];
    %     end
    %     filter_pars_2 = filter_pars;
    
    %=========================================================================
    %% MERGE SPOTS (if asked) that are to close too each other
    % spots that are closer than R_max will be merged into one spot and their data will be averaged
    
    %if spot_merge
%        [spots_clusters_all, cluster2_all] = cluster_Spots(spotData_1, R_merge);
 %       combineList = [cluster2_all(:,1), cluster2_all(:,7)];
 %       spotData_1 = merge_Spots(spotData_1, combineList);
  %  end
    
    %=========================================================================
    %% SHOW RESULTS of spot detection/filtering
    
    % frame to show
    %frame=1;
    
    % get data for specified channel...
 %   image = max(Fluo_1,[],3);
  %  spotData = spotData_1;
    
    % % get dada for specified channel...
    % image=max(Fluo_2,[],3);
    % %image=Fluo_1;
    % spotData=spotData_2;
    
    % setup the image scaling
%    im_min = min(image(:));
%    im_max = max(image(:));
%    im_range = im_max-im_min;
  %  im_0 = im_min+im_sc(1)*im_range;
 %   im_1 = im_min+im_sc(2)*im_range;
    
  %  if show_fig
  %      figure;
  %      imshow(image,[im_0, im_1],'InitialMagnification',200,'Border','tight');
  %      hold on
  %      if show_pars && ff==1
            disp('showing spot parameters: cyan = SpotIntensity/10000 red = IntensityRatio')
   %     end
    %    for ii=1:length(spotData)
    %        pos=spotData{ii}.spotPosition;
    %        SI=spotData{ii}.spotIntensity;
    %        IR=spotData{ii}.intensityRatio;
    %        plot(pos(1),pos(2),'+g','LineWidth',2)
    %        text(pos(1)+2,pos(2)+2,num2str(ii),'Color','g','FontSize',9)
    %        if show_pars
    %            text(pos(1)+1,pos(2)-3,num2str(SI/10000,3),'Color','c','FontSize',8)
     %           text(pos(1)-4,pos(2)-3,num2str(IR,3),'Color','r','FontSize',8)
     %       end
            
      %  end
    %end
    
    %if save_fig
     %   savefig(gcf,[path_3, tag,' ', filename_1(1:end-7),' Detected Spots.fig'],'compact')
     %   %print([save_fig_path,name_add,' R_Hist'],'-depsc')
    %end  
    
end % end of filenames loop
disp('DONE with All files!!!')


%% SAVE SPOTS and detection params

if ssave
    
%    save ([path_2, filename_out], 'spotDetection_params_1', 'spotDetection_params_2', 'path_1', 'filenames', 'filter_pars_2', 'filter_pars_1', 'filter_spots_1', 'filter_spots_2','dx', 'dz', 'spotData_1_All', 'spotData_2_All', '-v7')
   save ([path_2, filename_out], 'spotDetection_params_1', 'path_1', 'filenames', 'filter_pars_1', 'filter_spots_1', 'dx', 'dz', 'spotData_1_All', 'spot_merge', 'R_merge', '-v7')
   
    disp(['Data saved in ', [path_3, filename_out]])
    
end



%% NEXT SECTION




