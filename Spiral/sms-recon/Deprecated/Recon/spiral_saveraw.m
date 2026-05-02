function status = spiral_saveraw(recoInfo)
%function status = spiralraw(recoInfo)
%   recoInfo from setup_recon.m
  



%keyboard
kx = reshape(recoInfo.kx,[recoInfo.nro recoInfo.nl]);
ky = reshape(recoInfo.ky,[recoInfo.nro recoInfo.nl]);
nro = recoInfo.nro;
nroa = recoInfo.nroa;


skp = 128;


fid = fopen(recoInfo.fname_data);
bytes_to_skip_hdr = fread(fid,1,'uint32');
fseek(fid,bytes_to_skip_hdr,'bof');

bytesPerImage = recoInfo.num_coils*(((skp/2) * 2) + ((2*recoInfo.nroa) * 4)); 

%%%checking file size
%%fsz = (nroa*2*4+128)*nl*nsl*ntp*num_echo;


mkdir('Recon')
cd('Recon')
tmpvol = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl);
dat = zeros(recoInfo.nroa*recoInfo.num_adc,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils); 
dat = zeros(recoInfo.nro,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils); 


for imageIndex = 1:(recoInfo.ntp)
   for shotIndex = 1:recoInfo.nl %% load each interleaf's data into the buffer
        for sliceIndex = 1:recoInfo.nsl
            for adcIndex = 1:recoInfo.num_adc
                for coilIndex = 1:recoInfo.num_coils %% read in each coil data
 
                    %% Read and parse interleaf points
                   sMDH = ice_read_mdh_va21(fid);
                    [raw,rcount] = fread(fid,2*nroa,'float32');

                    if (rcount*recoInfo.num_adc < 2*nro)
		         sprintf('Reached end of file before running out of indices \n')
                         return
                    end

                    %% Put it in the correct buffer position
                        tmp = raw(1:2:2*nroa)+i*raw(2:2:2*nroa);
                        start_ind = (adcIndex-1)*nroa + 1;
                        end_ind = min(adcIndex*nroa,nro);
                        num_ind = end_ind - start_ind +1;
                        tmp2 = col(tmp(1:num_ind)).*exp(-i*2*pi*ky(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVy).*exp(-i*2*pi*kx(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVx);
                        dat(start_ind:end_ind,shotIndex,sliceIndex,coilIndex) = tmp2; 
                   
         	end % end coilIndex
      	     end  % end adcIndex
   	end  % end sliceIndex
    end  % end shotIndex


     save(sprintf('dat_%05d',imageIndex),'dat')
 


end % end imageIndex

    
status = 1;

cd ..

fclose('all')



