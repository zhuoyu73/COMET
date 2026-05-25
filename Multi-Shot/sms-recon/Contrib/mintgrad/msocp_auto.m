function [x,z,list_gap,list_dev]=msocp_auto(f,A,b,N,Nu,Kauto,maxgap,maxiter)

% [x,z,list_gap,list_dev]=msocp_auto(f,A,b,N,Nu,Kauto,maxgap,maxiter)
%
% all Matlab version of socp_auto()
%   only absolute tolerance stopping is implemented
%   list_gap and list_dev do not include gap and deviation of initial point
%   no time statistics
% Miguel S. Lobo -- 95/96


Mvar=0;		% 0: big-M can only grow; 1: big-M can grow or decrease

max_iter_plane=20;		% max. number of plane search iterations
max_lambda2=1e-2;		% stopping criterion for plane search
div_alpha=2;			% dividing const. for plane search
min_alpha=div_alpha^(-20);	% allow max. of 20 iteration of line search

[m,n]=size(A);
L=length(N);
% build cone variable selectors
SU=[];
ST=[];
T=[];
for i=1:L,
 SU=[SU; [zeros(N(i),i-1) [ones(N(i)-1,1); 0] zeros(N(i),L-i)]];
 ST=[ST; [zeros(N(i),i-1) [zeros(N(i)-1,1); 1] zeros(N(i),L-i)]];
 T=[T; zeros(N(i)-1,1); 1];	% or: T=sum(ST')';
end
SF=ST-SU;

% PRIMAL FEASIBLE POINT
x=zeros(n,1);
u=A*x+b;			% point in feasible hyperplane
k=-ST'*u+sqrt(SU'*(u.^2)+1);	% vector with "infeasibility" of each cone
x1=max(1,max(k));		% to make all feas. (min 1 => x1 feas.)
x=[x; x1];			% x1 will be added to all t(i)
M2=Kauto*T'*u;			% bound on max(t)

% DUAL FEASIBLE POINT
%	z=pinv(A')*f;  % (use this if A is square)
z=A'\f;				% pick any point satisfying hyperplane constr.
k=-ST'*z+sqrt(SU'*(z.^2)+1);	% vector with "infeasibility" of each cone
w2=max(1,max(k));		% to make all feas. (min 1 => w2 feas.)
z=z+T*w2;			% s(i)+=w2;
M1=Kauto*T'*z;			% bound on max(s)
w1=M1-T'*z;			% to satisfy extended hyperplane constr.
z=[z; w1; w2];			% extend dual variable

% EXTENDED PROBLEM
csum=T'*A;		% sum(c(i)')
dsum=T'*b;		% sum(d(i))
f=[f; M1];		% extend cost
A=[A T; zeros(1,n) 1; -csum 0]; %-L	% t(i)+=x1  and  v1=x1
b=[b; 0; (M2-dsum)];			% v2=M2-sum(t(i))
N=[N(:); 1; 1];
SU=[SU zeros(m,2); zeros(2,L) zeros(2,2)];
ST=[ST zeros(m,2); zeros(2,L) eye(2)];
SF=ST-SU;
T0=[T; 0; 0];
T1=[T; 1; 1];
m=m+2;
n=n+1;
L=L+2;

w=2*L+Nu*sqrt(2*L);	% duality gap vs. deviation from centrality weight

list_gap=[];
list_dev=[];
% list_zerr=[];

u=A*x+b;
gap=u'*z;
%	disp('extension done'); keyboard;
iter_out=0;
while (gap>maxgap) & (iter_out<maxiter),	% MAIN LOOP
 iter_out=iter_out+1;

 ff = -2*SF*(1./(SF'*(u.^2)));	% concatenation of
				%   -2/(t.^2-u'*u) * [-1;-1;...;-1;1]
				%   for each constraint
 gup = ff.*u;			% gradient of primal barrier wrt. u
 gu = gup + (w/gap)*z;		% gradient of potential wrt. u
 Hu = gup*gup'.*abs(SF*SF') + diag(ff);	% Hessian of primal barrier wrt. u
 Hx = A'*Hu*A;
 dx = - Hx\(A'*gu);		% primal search direction
%	dx = -(Hu*A)\gu;
 du = A*dx;
 dz = -(gu+Hu*du);		% residual provides  dual search direction
%	disp('directions found'); keyboard;

%	ddz = -A'\(A'*dz);
%	error = norm(ddz)/norm(dz)
%	dz = dz + ddz;

 % constants for plane search
 c1=gap;	c2=du'*z;	c3=u'*dz;
 d1=SF'*(u.^2);	d2=SF'*(u.*du);	d3=SF'*(du.^2);
 e1=SF'*(z.^2);	e2=SF'*(z.*dz);	e3=SF'*(dz.^2);
 p=0; q=0;		% initial value for search
 lambda2=inf;
 iter_plane=0;
 while (lambda2>max_lambda2) & (iter_plane<max_iter_plane),	% PLANE LOOP
  iter_plane=iter_plane+1;

  gp=w*c2/(c1+p*c2+q*c3)-2*sum((d2+p*d3)./(d1+p*2*d2+p^2*d3));
  gq=w*c3/(c1+p*c2+q*c3)-2*sum((e2+q*e3)./(e1+q*2*e2+q^2*e3));
  hp=-2*sum((d3.*(d1+p*2*d2+p^2*d3)-2*(d2+p*d3).^2)./(d1+p*2*d2+p^2*d3).^2);
  hq=-2*sum((e3.*(e1+q*2*e2+q^2*e3)-2*(e2+q*e3).^2)./(e1+q*2*e2+q^2*e3).^2);
  dp = -gp/hp;
  dq = -gq/hq;

  alpha=1;	p0=p;	q0=q;
  while 1,				% LINE SEARCH LOOP
   p=p0+alpha*dp;
   q=q0+alpha*dq;

   ua=u+p*du;
   za=z+q*dz;
   fi=[SF'*(ua.^2); ST'*ua; SF'*(za.^2); ST'*za];

   gpa=w*c2/(c1+p*c2+q*c3)-2*sum((d2+p*d3)./(d1+p*2*d2+p^2*d3));
   gqa=w*c3/(c1+p*c2+q*c3)-2*sum((e2+q*e3)./(e1+q*2*e2+q^2*e3));
   galpha=dp*gpa+dq*gqa;

   if all(fi>0) & galpha<=0; break; end	% successful exit of line search
   alpha=alpha/div_alpha;		% try smaller step
   if alpha<min_alpha			% too small - line search failure
    alpha=0;
    p=p0;
    q=q0;
    break;
   end
  end					% END OF LINE SEARCH
  lambda2 = hp*dp^2 + hq*dq^2;
%  lambda2 = hp*(alpha*dp)^2 + hq*(alpha*dq)^2;
 end					% END OF PLANE SEARCH LOOP

 x=x+p*dx;				% update primal and dual variables
 u=A*x+b;
 z=z+q*dz;

 M1new=Kauto*T0'*z;
 M2new=Kauto*(T0'*u-(L-2)*u(m-1));
 if (M1new>M1) | Mvar,
  M1=M1new;
  f(n)=M1;
  z(m-1)=(Kauto-1)*M1/Kauto;
 end;
 if (M2new>M2) | Mvar,
  M2=M2new;
  b(m)=M2-dsum;
  u(m)=(Kauto-1)*M2/Kauto;
 end;

 gap = u'*z;
 dev = -sum([log(SF'*(u.^2));log(SF'*(z.^2))]) + 2*L*(log(gap)-log(L));
% zerr = norm(A'*z-f) / norm(A'*z);

 list_gap=[list_gap; gap];
 list_dev=[list_dev; dev];
% list_zerr=[list_zerr; zerr];
end					% END OF MAIN LOOP

	disp('[alpha, beta]=');
	disp([u(m-1),z(m)]);
%	keyboard;

x=x(1:n-1);		% remove extra variables from solution
z=z(1:m-2);
