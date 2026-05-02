function vo = prepData(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y

if a.is.empty
	error empty
end
AA = a.A{1};
ntp = a.ntp;
V = AA.V;
num_coils = AA.num_coils; %size(sen,2);  % NOTE THIS IS RANK oF SVD Combination
rank_svd = AA.rank_svd;
vo = zeros(length(vi)*rank_svd/num_coils,1);
vi_index = 1;
vo_index = 1;
list=a.list_phase_encode;

for aa = 1:ntp
        kk=list(aa);
        sizeA = size(a.A{kk});
        vi_step = sizeA(1)*num_coils/rank_svd;
        vi_tmp = vi(vi_index:(vi_index+vi_step-1));
        vi_tmp = reshape(vi_tmp,[],num_coils);
        vo_tmp = vi_tmp*V(:,1:rank_svd);
        vo_step = numel(vo_tmp);
        vo(vo_index:(vo_index+vo_step-1)) = vo_tmp(:);
        vi_index = vi_index+vi_step;
        vo_index = vo_index +vo_step;
end

vo = vo(:);
