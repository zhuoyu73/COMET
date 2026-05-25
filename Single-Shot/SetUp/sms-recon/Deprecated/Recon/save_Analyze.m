function status = save_Analyze(str_stem, int_flag)
% function status = save_Analyze(str_stem, int_flag)
% str_stem = '' for all
% str_stem = 'sen' for only sen recons
% int_flag = 1 for interleaves

if ~exist('str_stem','var')
   str_stem = '';
end
if ~exist('int_flag','var')
   int_flag = 0;
end




load recoInfo
cd Recon
fnames = dir(sprintf('%s*.mat',str_stem));

for ii = 1:length(fnames)
    eval(sprintf('load %s',fnames(ii).name))
    if ~exist('img','var')
        eval(sprintf('img = %s;',str_stem))
    end
    name_to_use = fnames(ii).name;
    name_to_use = name_to_use(1:end-4);
    if int_flag
        [nx,ny,nz] = size(img);
        imgout = zeros(size(img));
        imgout(:,:,1:2:end) = img(:,:,1:round(nz/2));
        imgout(:,:,2:2:end) = img(:,:,round(nz/2)+1:end);
        img = imgout;
    end
    if isreal(img)
     write_ana(name_to_use,permute(img,[2 1 3]),recoInfo.szx, recoInfo.szy, recoInfo.szz,'FLOAT',1);
    else
     write_ana(name_to_use,permute(abs(img),[2 1 3]),recoInfo.szx, recoInfo.szy, recoInfo.szz,'FLOAT',1);
    end
end

cd ..
