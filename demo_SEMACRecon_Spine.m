clear
clc

%% Add path
addpath(genpath('./functions'));
addpath(genpath('./mrdCode'))

computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '/server/home/bli/MetalImaging/Code/SEMAC';
    ismrmrd_directory = '/server/home/bli/ismrmrd';
    addpath(src_directory);
elseif strcmp(computer_type,'MACA64')
    ismrmrd_directory = '/Users/bli/ismrmrd';
end
addpath(genpath(ismrmrd_directory));

%% Convert to ISMRMD 
data_folder = '/server/sdata/mri_data/disc/phantom/phantom0001_20231108/raw/SEMAC104/';
h5_SEMAC = [data_folder,'h5/meas_MID01939_FID12644_t2_tse_sag_semac8_p2_nofilter_104.h5'];
noise_SEMAC =  [data_folder,'noise/noise_meas_MID01939_FID12644_t2_tse_sag_semac8_p2_nofilter_104.h5'];

Mat_folder = [data_folder,'MAT/']; % Direcoty where MAT saved
kfill_savepath = fullfile(Mat_folder,'meas_MID01939_FID12644_t2_tse_sag_semac8_p2_nofilter_104_kfill.mat');
if ~exist(Mat_folder)
    mkdir(Mat_folder)
end


if ~exist(h5_SEMAC) || ~exist(noise_SEMAC)
        demo_convert_twix_to_ismrmrd([data_folder]);
end

%% Initialize processing flags
Read_flags.h5_fileList       = h5_SEMAC;
Read_flags.noise_fileList    = noise_SEMAC;
Read_flags.RemoveOS          = true; % remove oversampling
Read_flags.IgnoreSeg         = true; % concatanate segmemtns
Read_flags.DoAverage         = true; % Do averages (if 'average' was set during data acquistion)
Read_flags.CropPhaseEncoding = false;
Read_flags.Squeeze           = true;
Read_flags.os                = 2; % oversampling rate (Siemens Default value, don't change it)
Read_flags.noisePreWhitening = true;
Read_flags.GRAPPA            = true;
Read_flags.SliceOrder        = 'int';
Recon_flags.CoilComb         = 'grappa'; % 'grappa', ' walsh', 'sos'

%% Reading data
[c_img, kdata, noise, info] = readSEMAC_ismrmd(Read_flags); % kdata:[Nkx,Nky,Nkz,NSlice,NCoil]
[Nx, Ny, Nz, Ns, Ncoil] = size(kdata{1});
kdata = kdata{1};

%% Coil combination
zero_shift_idx = Nz/2; 
if Read_flags.GRAPPA || strcmp(Recon_flags.CoilComb, 'grappa') % GRAPPA to fill missing k-space
    grappa_start = tic;
    fprintf('Starting GRAPPA Recon:\n')
    kdata_fill = zeros(Nx,Ny,Nz,Ncoil,Ns);
    kernal_size = [5,5]; % kernal size
    parfor ns = 1:Ns
        kdata_2dcoils = squeeze(kdata(:,:,:,ns,:));
        kcalib = squeeze(kdata(:, info{1}.calib_idx, zero_shift_idx, ns, :));
        kdata_fill(:, :, :, :, ns) = GRAPPA_semsi(kdata_2dcoils, kcalib, kernal_size, 0.01);
        fprintf('Finish %d/%d slice: %.2f seconds \n', ns, Ns, toc(grappa_start))
    end
    fprintf('Total time: %2.f seconds', toc(grappa_start))
    kdata_fill = permute(kdata_fill,[1 2 3 5 4]); % [Nx,Ny,Nz,Ns,Ncoil]
    save(kfill_savepath, 'kdata_fill', '-v7.3');
    
    [Nx,Ny,Nz,Ns,Ncoil] = size(kdata_fill);

    % Walsh coil estimation following GRAPPA:
    csm_start = tic;
    fprintf('Estimating coil sensitvity...')
    csm_option.method = 'walsh'; % coil estimation method: 'walsh', 'sos'
    csm_option.cal_shape = [24 24]; %calibration region
    csm_option.kdata = squeeze(kdata_fill(:,:,zero_shift_idx,:,:));
    [csm, cal_im] = coil_estimation(csm_option);
    csm = repmat(csm,[1 1 1 1 Nz]); % [Nx,Ny,Ns,Ncoil,Ne]
    csm = permute(csm,[1 2 5 3 4]); % [Nx,Ny,Ne,Ns,Ncoil]
    SEMAC_coils = fft3c(kdata_fill); % [Nx,Ny,Ne,Ns,Ncoil]
    SEMAC_CC =  sum( conj(csm) .* SEMAC_coils, 5)  ./ sqrt(sum(abs(csm).^2,5));
    fprintf('Done!!\b %f s \n',toc(csm_start))
    
    
    
elseif strcmp(Recon_flags.CoilComb, 'walsh') % Adpative coil combination
    SEMAC_coils = fft3c(kdata); 
    csm_option.method = 'walsh'; % coil estimation method: 'walsh', 'sos'
    csm_option.cal_shape = [32 32]; %calibration region
    csm_option.kdata = squeeze(kdata{1}(:,:,4,:,:)); % Please change this number to the kz of 0. (7 for SEMAC12; 4 for SEMAC6)

    
    [csm, cal_im] = coil_estimation(csm_option);
    
    csm = repmat(csm,[1 1 1 1 Nz]);
    csm = permute(csm,[1 2 5 3 4]);
    
    SEMAC_CC =  sum( conj(csm) .* SEMAC_coils, 5)  ./ sqrt(sum(abs(csm).^2,5)); % Coil combined images ./ sqrt(sum(abs(csm).^2,4))

elseif strcmp(Recon_flags.CoilComb, 'sos')
    SEMAC_CC = sos(fft3c(kdata));
end
return

%% Crop PE
% dset = ismrmrd.Dataset(h5_SEMAC, 'dataset');
% header = ismrmrd.xml.deserialize(dset.readxml);
% pe_range = header.encoding.reconSpace.matrixSize.y/2;
%  if  header.encoding.encodingLimits.kspace_encoding_step_1.maximum >= pe_range
%             remain_start = header.encoding.encodingLimits.kspace_encoding_step_1.center - pe_range/2 + 1;
%             remain_end   = header.encoding.encodingLimits.kspace_encoding_step_1.center + (pe_range - 1 -pe_range/2) + 1;
%             kdata_crop = kdata_fill(:,remain_start:remain_end,:,:,:);
%             SEMAC_CC_Crop = sos(fft3c(kdata_crop));
%  end
%         
%% Re-order slice
SEMAC_CC_correct_slice = reorder_slice(SEMAC_CC, Read_flags, Ns);

%% SEMAC  combination
kz_s_idx = repmat([(1:Nz).';zeros(Ns,1)],[1 Ns]);
for ns = 1:Ns
    kzs_map(:,ns) = circshift(kz_s_idx(:,ns),ns-1);
end
kzs_map = flip(kzs_map,2);
center_z = ceil(Nz/2)+1;
image_flip = flip(SEMAC_CC_correct_slice,4);
SEMAC_Recon_sos_2 = zeros(size(SEMAC_CC_correct_slice,1,2,4));

for ss = 1:Ns
    centerZ_in_map = find(kzs_map(:,ss) == center_z);
    slice_idx = find(kzs_map(centerZ_in_map,:));
    z_idx     = kzs_map(centerZ_in_map,slice_idx);
    for nsum = 1:length(z_idx)
        SEMAC_Recon_sos_2(:,:,ss)  = SEMAC_Recon_sos_2(:,:,ss) + abs(image_flip(:,:,z_idx(nsum),slice_idx(nsum))).^2;
    end
end
SEMAC_Recon_sos = sqrt(SEMAC_Recon_sos_2);

%% SEMAC combination (SEMAC=6;)
SEMAC_Recon_sos = zeros(size(SEMAC_CC_correct_slice,1,2,4));
for ss = 1:size(SEMAC_Recon_sos,3)
    if (ss-2 >0) & (ss+3 <= Ns) 
    SEMAC_Recon_sos(:,:,ss) = (sqrt(abs(SEMAC_CC_correct_slice(:,:,6,ss-2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,5,ss-1)).^2 +...
        abs(SEMAC_CC_correct_slice(:,:,4,ss)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,3,ss+1)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,2,ss+2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,1,ss+3)).^2));
    end
end

SEMAC_Recon_complex = zeros(size(SEMAC_CC_correct_slice,1,2,4));
for ss = 1:size(SEMAC_Recon_complex,3)
    if (ss-2 >0) & (ss+3 <= Ns) 
    SEMAC_Recon_complex(:,:,ss) = SEMAC_CC_correct_slice(:,:,6,ss-2)+...
        SEMAC_CC_correct_slice(:,:,5,ss-1) +...
        SEMAC_CC_correct_slice(:,:,4,ss)+...
        SEMAC_CC_correct_slice(:,:,3,ss+1)+...
        SEMAC_CC_correct_slice(:,:,2,ss+2)+...
        SEMAC_CC_correct_slice(:,:,1,ss+3);
    end
end

return
%% SEMAC combination (SEMAC=8)
SEMAC_Recon_sos = zeros(size(SEMAC_CC_correct_slice,1,2,4));
for ss = 1:size(SEMAC_Recon_sos,3)
    if (ss-3 >0) && (ss+4 <= Ns) 
    SEMAC_Recon_sos(:,:,ss) = sqrt(abs(SEMAC_CC_correct_slice(:,:,8,ss-3)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,7,ss-2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,6,ss-1)).^2 +...
        abs(SEMAC_CC_correct_slice(:,:,5,ss)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,4,ss+1)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,3,ss+2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,2,ss+3)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,1,ss+4)).^2);
    end
end
%% Transpose 
SEMAC_Recon_sos_rot = zeros(size(SEMAC_Recon_sos,2,1,4));
for ss = 1:Ns
    SEMAC_Recon_sos_rot(:,:,ss) = SEMAC_Recon_sos(:,:,ss).';    
end

return

%%
ss12_kz = cat(3,(SEMAC_CC_correct_slice(:,:,8,ss-3)),...
        (SEMAC_CC_correct_slice(:,:,7,ss-2)),...
        (SEMAC_CC_correct_slice(:,:,6,ss-1)),...
        (SEMAC_CC_correct_slice(:,:,5,ss)),...
        (SEMAC_CC_correct_slice(:,:,4,ss+1)),...
        (SEMAC_CC_correct_slice(:,:,3,ss+2)),...
        (SEMAC_CC_correct_slice(:,:,2,ss+3)),...
        (SEMAC_CC_correct_slice(:,:,1,ss+4)));
    



%% SEMAC combination (SEMAC=12 ; 21 slices)
SEMAC_Recon_sos = zeros(size(SEMAC_CC,1,2,4));
for ss = 1:21
    if (ss-5 >0) && (ss+ 6<= 21) 
    SEMAC_Recon_sos(:,:,ss) = sqrt(abs(SEMAC_CC(:,:,12,ss-5)).^2+...
        abs(SEMAC_CC(:,:,11,ss-4)).^2+...
        abs(SEMAC_CC(:,:,10,ss-3)).^2+...
        abs(SEMAC_CC(:,:,9,ss-2)).^2+...
        abs(SEMAC_CC(:,:,8,ss-1)).^2 +...
        abs(SEMAC_CC(:,:,7,ss)).^2+...
        abs(SEMAC_CC(:,:,6,ss+1)).^2+...
        abs(SEMAC_CC(:,:,5,ss+2)).^2+...
        abs(SEMAC_CC(:,:,4,ss+3)).^2+...
        abs(SEMAC_CC(:,:,3,ss+4)).^2+...
        abs(SEMAC_CC(:,:,2,ss+5)).^2+...
        abs(SEMAC_CC(:,:,1,ss+6)).^2);
    end
end


%% Save results
% 1) Save full k-space
pattern = 'h5/(.*?)\.h5';
matches = regexp(h5_SEMAC, pattern, 'tokens')
output_fullkspace = [matches{1}{1},'_Fullkspace'];
if  Read_flags.GRAPPA || strcmp(Recon_flags.CoilComb, 'grappa')
    time_save = tic;
    fprintf('Saving Full k-space in %s\n',[data_folder, output_fullkspace])
    save([data_folder, output_fullkspace],'kdata_full','-v7.3');
    fprintf('Done! %.2f seconds \n', toc(time_save))
end 
% 2) Save Recon images
output_recon_img = [matches{1}{1},'img_semac6_p2_grappspe']
fprintf('Saving SEMSI images in %s\n',[data_folder,output_recon_img])
save([data_folder,output_recon_img],'SEMAC_Recon', '-v7.3');
fprintf('Done! %.2f seconds \n', toc(time_save))