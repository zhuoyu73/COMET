 function vo = mtimes(a, vi)
%	MRI "forward projection" y=A*x and backprojection x = (A')*y
%   Brad Sutton, Univ. Michigan, June 2002

if a.is.empty
	error empty
end

v = a.v;
R = a.R;
ntp = a.ntp;
sizeA = size(a);
if ~a.is.transpose 
    vi = reshape(vi,[],R);
    vi_tmp = vi*v;
    vo = zeros(sizeA(1),1);  
    start_index = 1;
    for aa = 1:ntp
        A_tmp = a.A{aa};
        sizeA = size(A_tmp);
        lengthA_tmp = sizeA(1);
        vo(start_index:(lengthA_tmp+start_index-1)) =  A_tmp*col(vi_tmp(:,aa)); 
        start_index = start_index+lengthA_tmp;
    end
    vo = vo(:);
else        
    vo_tmp = zeros(sizeA(2)/R, ntp);
    start_index = 1;
    for aa = 1:ntp
        A_tmp = a.A{aa};
        sizeA = size(A_tmp);
        lengthA_tmp = sizeA(1);
        vi_tmp = vi(start_index:(lengthA_tmp+start_index-1));  
        start_index = start_index+lengthA_tmp;
        vo_tmp(:,aa) = A_tmp'*vi_tmp;        
    end
    vo = vo_tmp*v';
    vo = vo(:);
end





