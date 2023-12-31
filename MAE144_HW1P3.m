clear; close all; clc;

ys=[20 20]; xs=[1 10]; h=0.01; Ds=RR_tf(ys,xs); omegac=sqrt(10); [Dz]=DL_C2D_matched(Ds,h,omegac)

%% MAE144 HW1 Problem3
function [dz] = DL_C2D_matched(ds,h, omgBar) 
% Matched pole-zero method that transform continuous time transfer
% function to descrete time transfer function.
% Ds - Continuouse time transfer function
% omgBar - frequency of interest


    if ~isa(h,'RR_poly'), h = 0.001; end % estimate a small value time step
    if ~isa(omgBar,'RR_poly'), omgBar = 0; end % check if omega bar input exists
    % h = 0.001;
    % omgBar = 0;
    
    zeros = RR_roots(ds.num); % finding the zeros and poles of continuous tf
    poles = RR_roots(ds.den);
    z_zeros = exp(zeros .* h); % initial s-to-z conversion
    z_poles = exp(poles .* h);

    m = ds.num.n; % check infinite zeros' existance
    n = ds.den.n;
    infz = m-n;
    if infz > 0 % including infinite zeros into the numerator, but include one less to make tf proper
        for j = 1:(infz-1), z_zeros = [z_zeros -1]; end
    end

    prodnum = 1; % product of numerator
    for j = 1:numel(z_zeros) prodnum = prodnum * (h * omgBar + z_zeros(j)); end
    prodden = 1; % product of denominator
    for j = 1:numel(z_poles) prodden = prodden * (h * omgBar + z_poles(j)); end
    kfac = prodnum/prodden; % finding the correct gain

    dz = kfac * RR_tf(z_zeros,z_poles); % return the fully transformed descrete transfer function
    dz.h = h;

end


% function [Dz]=RR_C2D_tustin(Ds,h,omegac)
% % function [Dz]=RR_C2D_tustin(Ds,h,omegac)
% % Convert Ds(s) to Dz(z) using Tustin's method.  If omegac is specified, prewarping is applied
% % such that the dynamics of Ds(s) in the vicinity of this critical frequency are mapped correctly.
% % TEST: ys=20*[1 1]; xs=[1 10]; h=0.01; Ds=RR_tf(ys,xs); omegac=sqrt(10); [Dz]=RR_C2D_tustin(Ds,h,omegac)
% %       disp('Corresponding Matlab solution:')
% %       opt = c2dOptions('Method','tustin','PrewarpFrequency',omegac); c2d(tf(ys,xs),h,opt)
%     if nargin==2, f=1; else, f=2*(1-cos(omegac*h))/(omegac*h*sin(omegac*h)); end
%     c=2/(f*h); m=Ds.num.n; n=Ds.den.n; b=RR_poly(0); a=b;
%     fac1=RR_poly([1 1]); fac2=RR_poly([1 -1]);
%     for j=0:m; b=b+Ds.num.poly(m+1-j)*c^j*fac1^(n-j)*fac2^j; end
%     for j=0:n; a=a+Ds.den.poly(n+1-j)*c^j*fac1^(n-j)*fac2^j; end, Dz=RR_tf(b,a); Dz.h=h;
% end % function RR_C2D_tustin

