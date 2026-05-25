function [tformList] = calc2Dregistrations_ZS(imgNav)
% Authors:
% Alex Cerjanic & Zhuoyu Shi
% June 2025

sizes = size(squeeze(imgNav));
Nx = sizes(1);
Ny = sizes(2);
Nz = sizes(3);
Nshots = sizes(4);
Npartitions = sizes(5);
Nslabs = sizes(6);
Nimages = sizes(7);

navImages = squeeze(abs(imgNav));

tformList = zeros(3,3,Nshots,Npartitions,Nslabs,Nimages);

sliceToUse = 3;

for image = 1:Nimages
    for slab = 1:Nslabs
        refImage = squeeze(navImages(:,:,sliceToUse,1,1,slab,1));
        for par = 1:Npartitions
            for shot = 1:Nshots
                imgToRegister = squeeze(abs(navImages(:,:,sliceToUse,shot,par,slab,image)));
                tform = imregcorr(imgToRegister,refImage);
                tformList(:,:,shot,par,slab,image) = tform.A;
                tformList(1:2,3,shot,par,slab,image) = tformList(1:2,3,shot,par,slab,image)*3;
            end
        end
    end
end

tformList = reshape(tformList,3,3,1,2,2,16,1,1,1,24);
tformList = repmat(tformList,1,1,4,1,1,1,1,1,1,1);
affine_matrices = tformList;
save(fullfile(pwd,'affine_matrices_crosscor_ZS.mat'),'affine_matrices')
end
