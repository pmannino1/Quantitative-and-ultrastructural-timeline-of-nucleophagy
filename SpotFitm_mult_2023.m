function [All_g_spotData, master_double] = SpotFitm_mult_2023(spots_path, images_path, output_path)

% How to run function:
% 1. Load all three variables (images_path = '' etc., one-by-one in command line)
%    Then, in command line:


   % [All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cell 1 spot data.xlsx', '\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Images\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cell 1 gaussians.xlsx')

    %[All_g_spotData, master_double] = SpotFitm_mult_2023('\Users\pm695\Desktop\Lusk Lab\Data\GFP counting experiments\05242023\FFC corrected images\Output\WT\WT Atg39-GFP spot data_rep 2.xlsx', '\Users\pm695\Desktop\Lusk Lab\Data\GFP counting experiments\05242023\FFC corrected images\WT\', '\Users\pm695\Desktop\Lusk Lab\Data\GFP counting experiments\05242023\FFC corrected images\Output\WT\Atg39-GFP gaussians WT_rep2\Cell 1 gaussians.xlsx')


   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cell 1 spot data.xlsx', '\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Images\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cell 1 gaussians.xlsx')

   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cropped output\Cell 1 spot datay_cropped.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Images\Images_2\cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cropped output\Cell 1 gaussians_cropped.xlsx')

   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cropped output\Cell 1 cropped spot data for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Images\Images_2\cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\0\Output\Cropped output\Cell 1 gaussians_cropped_2.xlsx')

   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\output\Focus 1\Cell 3 spot data for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\Images\Focus 1\Cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\output\Focus 1\Cell 3 gaussians_croppedfocus1.xlsx')

   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\output\Focus 2\Cell 3 spot data_for gaussians focus2.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\Images\Focus 2\Cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\2\output\Focus 2\Cell 3 gaussians_croppedfocus2.xlsx')

   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\3\output\second event\Cell 4 spot data_event2_for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\3\Images\Images 2\Second event\Cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\3\output\second event\Cell 4 gaussians_event2.xlsx')
   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\4\output\Cell 5 spot data_for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\4\Images\Images 2\cropped\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06162023\Output ROIs\01\4\output\Cell 5 gaussians.xlsx')

      %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\03272023\ROIs\02\12\output\Cell 12 spot data_later_for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\03272023\ROIs\02\Later\112\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\03272023\ROIs\02\12\output\Cell 12_Gaussians_later.xlsx')



   %[All_g_spotData, master_double] = SpotFitm_mult_2023('C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06172023\ROIs\01\6\output\Cell 6 spot data_later_for gaussians.xlsx', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06172023\ROIs\01\Later\6\', 'C:\Users\pm695\Desktop\Lusk Lab\Data\LLSM\06172023\ROIs\01\6\output\Cell 6 Gaussians_later.xlsx')


% INPUTS:
% spots_path: A PATH to a csv file containing the particles to fit
% Gaussians & G_ERFs to. One image, indexed in the seventh row, may have
% multiple spots
% images_path: A PATH to a folder containing *only* the .tif files to be
% considered. 
% output_path: A PATH to a proposed file (with extension) where spot data 
% will be exported. 

% OUTPUTS:
% all_g_spotData: a 1 x tracks CELL ARRAY with each cell containing the
% output of SpotDetection_2023 or findIrregularSpots. Each cell in all_g_SD
% contains a 1xframes cell array corresponding to one track. Each cell in
% that array contains a 1x1 struct containing a particular frame's
% position, intensity, intensity ratio, frame #, Fluo_1, and results of
% fitSpots_2023.m fitting + parameters.
% master_double: written to the output path. Columns: x-coord, y-coord,
% z-coord, z-coord, frame, volume of Gaussian fit, spot intensity,
% intensity ratio

% Required Functions:
% Get_Filenames
% fitSpots_2023
% load_tiff_stack

spots = readtable(spots_path);
G = findgroups(spots{:,7});
Tc = splitapply( @(varargin) varargin, spots, G);

subTables = cell(size(Tc, 1));
% Create sub tables
for i = 1:size(Tc, 1)
    subTables{i} = table(Tc{i, :}, 'VariableNames', spots.Properties.VariableNames);
end

All_g_spotData = {};

for st=1:size(subTables)

    singleframe_spotData = {};

    for rows=1:height(subTables{st}) % go through single track rows

        xyzPosition = table2array(subTables{st}(rows,1:3));
        intRat = table2array(subTables{st}(rows,6));
        spotInt = table2array(subTables{st}(rows,5));
        spotFrame = table2array(subTables{st}(rows,7));

        % make struct for each row
        particle.spotPosition = xyzPosition;
        particle.intensityRatio = intRat;
        particle.spotIntensity = spotInt; 
        particle.csvFrame = spotFrame;

        % assign struct to a cell, input to singletrack data
        singleframe_spotData{rows} = {particle};

    end 

    All_g_spotData{st} = singleframe_spotData;

end

%% Loading Fluo_1 and Gaussain results to each particle in All_g_spotData

for oneframe=1:length(All_g_spotData)

    % load image files in order
    tif_filenames = Get_Filenames(images_path);

    % load the single frame
    [Fluo_1, image_info]=load_tiff_stack(strcat(images_path,tif_filenames{oneframe}));

    for particles=1:width(All_g_spotData{oneframe})
        particlestruct = All_g_spotData{oneframe}{particles}{1};
        particlestruct.Fluo_1 = Fluo_1;

        % fitting + assigning
        [params, FitResults_G, FitResults_GE] = fitSpots_2023(particlestruct.Fluo_1, particlestruct);
        particlestruct.FitResults_G = FitResults_G;
        particlestruct.FitResults_GE = FitResults_GE;
        particlestruct.params = params;

        b = particlestruct.FitResults_G.coeff(4);
        sigma2 = particlestruct.FitResults_G.coeff(5);
        c = particlestruct.FitResults_G.coeff(1);

        Volume = (2 * b * sigma2 * pi);

        particlestruct.gaussVol = Volume;

        All_g_spotData{oneframe}{particles}{1} = particlestruct;

    end

end

%% Compiling a double of each Gaussian volume, SI, IR for the dataset

master_double = [];

for frame=1:length(All_g_spotData)
    zeroframe = zeros(width(All_g_spotData{frame}), 7);
    for rows=1:height(zeroframe)
        particlestruct = All_g_spotData{frame}{rows}{1};
        zeroframe(rows, 1) = particlestruct.spotPosition(1, 1);
        zeroframe(rows, 2) = particlestruct.spotPosition(1, 2);
        zeroframe(rows, 3) = particlestruct.spotPosition(1, 3);
        zeroframe(rows, 4) = particlestruct.csvFrame;
        zeroframe(rows, 5) = particlestruct.gaussVol;
        zeroframe(rows, 6) = particlestruct.spotIntensity;
        zeroframe(rows, 7) = particlestruct.intensityRatio;
    end

    master_double = cat(1,master_double,zeroframe);
    
end

writematrix(master_double, output_path);
disp("Particle Table exported.")

end