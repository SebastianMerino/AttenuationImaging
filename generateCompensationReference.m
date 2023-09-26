%%
clear,clc
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\', ...
    'Attenuation\DataQUS_4_Merino'];
refDir = [baseDir,'\References\P4-CUELLO-2'];
refFiles = dir([refDir,'\*.mat']);
load([refDir,'\',refFiles(1).name])

%% Diffraction compensation
att_ref_dB = [0 0.3 0]; % Prior attenuation
att_ref = (att_ref_dB(1)*(f.^2) + att_ref_dB(2)*f + att_ref_dB(3) )/8.686;

% Attenuation reference
att_ref_map = zeros(m,n,r);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Four 'samples' of the reference phantom
Sp_ref = zeros(m,n,r);
Sd_ref = zeros(m,n,r);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        % Reference 1
        sub_block_p = ref1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp1,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd1,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 2
        sub_block_p = ref2(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref2(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp2,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd2,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 3
        sub_block_p = ref3(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref3(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp3,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd3,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 4
        sub_block_p = ref4(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref4(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp4,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd4,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);

        tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
        tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));
    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;
