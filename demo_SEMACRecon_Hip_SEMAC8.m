clear
clc

%% Add path
addpath(genpath('./functions'));
addpath(genpath('./mrdCode'))

computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '/server/home/bli/MetalImaging/Code/SEMSI';
    ismrmrd_directory = '/server/home/bli/ismrmrd';
    data_parent_directory = '';
elseif strcmp(computer_type,'MACI64')
    ismrmrd_directory = '/Users/bli/ismrmrd';
end
addpath(genpath(ismrmrd_directory));
addpath(genpath(src_directory));

%% Convert to ISMRMD 
data_folder = '/server/home/bli/Data/SEMSI/20230927_SEMIS_vol714_RightHipTi/SEMAC/';
h5_SEMAC = [data_folder,'h5/meas_MID00459_FID08154_pd_tse_cor_semac8_p2.h5'];
noise_SEMAC =  [data_folder,'noise/noise_meas_MID00459_FID08154_pd_tse_cor_semac8_p2.h5'];

if ~exist(h5_SEMAC) || ~exist(noise_SEMAC)
        demo_convert_twix_to_ismrmrd([data_folder]);
end

%% Initialize processing flags
Read_flags.h5_fileList       = h5_SEMAC;
Read_flags.noise_fileList    = noise_SEMAC;
Read_flags.RemoveOS          = true; % remove oversampling
Read_flags.IgnoreSeg         = true; % concatanate segmemtns
Read_flags.DoAverage         = true; % Do averages (if 'average' was set during data acquistion)
Read_flags.CropPhaseEncoding = true;
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
if Read_flags.GRAPPA || strcmp(Recon_flags.CoilComb, 'grappa') % GRAPPA to fill missing k-space
    grappa_start = tic;
    fprintf('Starting GRAPPA Recon:\n')
    kdata_full = zeros(Nx,Ny,Nz,Ncoil,Ns);
    zero_sp_idx = 4; % zero spectral dimension index
    kernal_size = [5,5]; % kernal size
    parfor ns = 1:Ns
        kdata_2dcoils = squeeze(kdata(:,:,:,ns,:));
        kcalib = squeeze(kdata(:, info{1}.calib_idx, zero_sp_idx, ns, :));
        kdata_full(:, :, :, :, ns) = GRAPPA_spe(kdata_2dcoils, kcalib, kernal_size, 0.01);
        fprintf('Finish %d/%d slice: %.2f seconds \n', ns, Ns, toc(grappa_start))   
    end
    fprintf('Total time: %2.f seconds', toc(grappa_start))
    kdata_full = permute(kdata_full,[1 2 3 5 4]);
    SEMAC_coils = fft3c(kdata_full);
    SEMAC_CC = sos(SEMAC_coils);
    
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
dset = ismrmrd.Dataset(h5_SEMAC, 'dataset');
header = ismrmrd.xml.deserialize(dset.readxml);
pe_range = header.encoding.reconSpace.matrixSize.y/2;
 if  header.encoding.encodingLimits.kspace_encoding_step_1.maximum >= pe_range
            remain_start = header.encoding.encodingLimits.kspace_encoding_step_1.center - pe_range/2 + 1;
            remain_end   = header.encoding.encodingLimits.kspace_encoding_step_1.center + (pe_range - 1 -pe_range/2) + 1;
            kdata_crop = kdata_full(:,remain_start:remain_end,:,:,:);
            SEMAC_CC_Crop = sos(fft3c(kdata_crop));
 end
        
%% Re-order slice

SEMAC_CC_correct_slice = reorder_slice(SEMAC_CC, Read_flags, Ns);
SEMAC_CC_Crop_correct_slice = reorder_slice(SEMAC_CC(:,remain_start:remain_end,:,:,:), Read_flags, Ns);


%% SEMAC combination (SEMAC=8)
SEMAC_Recon = zeros(size(SEMAC_CC_Crop_correct_slice,1,2,4));
for ss = 1:size(SEMAC_Recon,3)
    if (ss-3 >0) && (ss+4 <= Ns) 
    SEMAC_Recon(:,:,ss) = (sqrt(abs(SEMAC_CC_Crop_correct_slice(:,:,8,ss-3)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,7,ss-2)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,6,ss-1)).^2 +...
        abs(SEMAC_CC_Crop_correct_slice(:,:,5,ss)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,4,ss+1)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,3,ss+2)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,2,ss+3)).^2+...
        abs(SEMAC_CC_Crop_correct_slice(:,:,1,ss+4)).^2));
    end
end

return

%%
ss8_kz = cat(3,abs(SEMAC_CC_correct_slice(:,:,8,ss-3)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,7,ss-2)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,6,ss-1)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,5,ss)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,4,ss+1)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,3,ss+2)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,2,ss+3)).^2,...
        abs(SEMAC_CC_correct_slice(:,:,1,ss+4)).^2);
    
%% SEMAC combination (SEMAC=6 ;15 slices)
SEMAC_Recon = zeros(size(SEMAC_CC_correct_slice,1,2,4));
for ss = 1:size(SEMAC_Recon,3)
    if (ss-2 >0) & (ss+3 <= Ns) 
    SEMAC_Recon(:,:,ss) = rot90(sqrt(abs(SEMAC_CC_correct_slice(:,:,6,ss-2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,5,ss-1)).^2 +...
        abs(SEMAC_CC_correct_slice(:,:,4,ss)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,3,ss+1)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,2,ss+2)).^2+...
        abs(SEMAC_CC_correct_slice(:,:,1,ss+3)).^2));
    end
end

return


%% SEMAC combination (SEMAC=12 ; 21 slices)
SEMAC_Recon = zeros(size(SEMAC_CC,1,2,4));
for ss = 1:21
    if (ss-5 >0) && (ss+ 6<= 21) 
    SEMAC_Recon(:,:,ss) = sqrt(abs(SEMAC_CC(:,:,12,ss-5)).^2+...
        abs(SEMAC_CC(:,:,11,ss-4)).^2+...
        abs(SEMAC_CC(:,:,10,ss-3)).^2+...
        abs(SEMAC_CC(:,:,9,ss-2)).^2+...
        abs(SEMAC_CC(:,:,8,ss-1)).^2 +...
        abs(SEMAC_CC(:,:,7,ss)).^2+...
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
output_recon_img = [matches{1}{1},'img_semac8_p2_grappa']
fprintf('Saving SEMSI images in %s\n',[data_folder,output_recon_img])
save([data_folder,output_recon_img],'SEMAC_Recon','SEMAC_CC_correct_slice','header', '-v7.3');
fprintf('Done! %.2f seconds \n', toc(time_save))