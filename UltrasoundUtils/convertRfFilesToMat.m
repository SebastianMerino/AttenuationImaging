function [] = convertRfFilesToMat(Directory)
% Converts all .rf files to .mat files in the directory

rfFiles = dir([Directory,'\*.rf']);
for iAcq = 1:length(rfFiles)
    out =lectura_OK([Directory,'\',rfFiles(iAcq).name]);
    save([Directory,'\',rfFiles(iAcq).name(1:end-3)],'-struct','out');
end

end