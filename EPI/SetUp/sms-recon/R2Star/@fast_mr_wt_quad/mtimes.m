 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
	error empty
end

G = a.G;
% epiflg = isa(G,'Gfft');
i = sqrt(-1);
tt = a.tt;
N=a.N;
L = size(a.int,1);
we = a.we;
Int = a.int;
TE = min(tt);
tt = tt-TE;
nz = size(we,3);



N1 = size(we,1);
N2 = size(we,2);
szz = size(a);
N3 = szz(1);
% N4 = szz(2);
% stx = (sqrt(N3)-sqrt(N2))/2;
% r2 = reshape(a.HH*a.r2(:)./nz,N1,N2); 
r2=a.r2;
if ~(size(tt) == size(a.int,2))
  sprintf('Number of time points not the same between P and Int')
  keyboard
end

if (L == 1)
   tau = 0;
else
  tau = (max(tt)-min(tt)+eps)/(L-1);
end



if ~a.is.transpose
	vi = a.HH*vi(:)./nz;   % oversampled in either direction
	vo = zeros(size(tt(:)));
    
    for ii = 1:nz
        wec = squeeze(we(:,:,ii));
        r2c = squeeze(r2(:,:,ii));


        vic = exp(-(i*wec +r2c)*TE).*reshape(vi(:),N1,N2);

        for ll = 1:L
            Wo = exp(-(i*wec+r2c).*((ll-1)*tau));
            aa = Int(ll,:).';

            vo = vo + aa.*(G*(Wo.*vic));

        end

    end
    
	if a.flgswt
		vo = vo.*a.swt;
	end
		
	    
else
    
	if a.flgswt
		vi = vi.*a.swt;  % Transpose of real sinc-weighting
    end
    
	vo = zeros(size(we,1),size(we,2));
    
	for ii = 1:nz
        
		wec = squeeze(we(:,:,ii));
		r2c = squeeze(r2(:,:,ii));
        
		for ll = 1:L
			Wo = exp((i*conj(wec)-r2c).*((ll-1)*tau));
			aa = Int(ll,:)';
            
			if ll == 1
				voc = Wo.*reshape(G'*(aa.*vi(:)),N1,N2);
			else
				voc = voc + Wo.*reshape(G'*(aa.*vi(:)),N1,N2);
            end
            
		end
		vo = vo +exp((i*conj(wec)-r2c).*TE).*voc;
        
	end
		vo = a.HH'*vo(:)./nz;
end





