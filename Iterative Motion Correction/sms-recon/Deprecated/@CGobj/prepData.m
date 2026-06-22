 function vo = prepData(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

if a.is.empty
	error empty
end

V = a.V;
num_coils = a.nc; %size(sen,2);  % NOTE THIS IS RANK oF SVD Combination
rank_svd = a.coil_rank;
nn = size(a.A{1});
nx = nn(1);
ny = nn(2);
nshots = a.n;

if (size(vi,1) ~= (num_coils*nx*nshots))
   sprintf('Error in prepData for coil combination SVD: sense_svd \n')
   return
end

     
vi = reshape(vi,nx, num_coils, nshots);

vo = zeros(nx,rank_svd,nshots);
for ll = 1:nshots
    vo(:,:,ll) = vi(:,:,ll)*V(:,1:rank_svd);
end

vo = vo(:);

