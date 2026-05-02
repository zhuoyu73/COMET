function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
	error empty
end

v = a.v;
R = a.R;
ntp = a.ntp;
list=a.A_index_list;
sizeA = size(a);
if ~a.is.transpose 
    vi = reshape(vi,[],R);
    vi_tmp = vi*v;
    vo = zeros(sizeA(1),1);  
    for aa = 1:ntp
        size_Atmp = size(A_tmp);
	    kk=list(aa);
        A_tmp = a.A{kk};
       vo((aa-1)*size_Atmp(1)+1:aa*size_Atmp(1)) =  A_tmp*col(vi_tmp(:,aa)); 
    end
    vo = vo(:);
else        
    vo_tmp = zeros(sizeA(2)/R, ntp);
    for aa = 1:ntp
        kk=list(aa);
        A_tmp = a.A{kk};
        size_Atmp = size(A_tmp);
        vi_tmp = vi((aa-1)*(size_Atmp(1))+1:aa*(size_Atmp(1)));  
        vo_tmp(:,aa) = A_tmp'*vi_tmp;        
    end
    vo = vo_tmp*v';
    vo = vo(:);
end

