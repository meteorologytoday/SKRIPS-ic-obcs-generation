function [To]=xyexpand(Tin,niter)
%% Tout=xyexpand(Tin,3);

%% This function takes 2D input and replaces any NaN
%% neighbouring any real numbers with the average of
%% those numbers. This procedure is repeated niter
%% times. Note that the ocean/land volume ratio expands
%% rapidly with each iteration.

%% In short, this is a gap-filler and extrapolator of
%% rudimentary order applied to horizontal slices.

%% Created  08/15/99 by adcroft@mit.edu
%% Modified 11/11/99 by adcroft@mit.edu
%% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

[nx,ny]=size(Tin);
ip=[2:nx 1];
im=[nx 1:nx-1];
jp=[2:ny ny];
jm=[1 1:ny-1];

T=Tin;
for iter=1:niter,
    
    msk=ones(nx,ny);
    L=isnan(T);
    msk( L )=0;
    
    To=T;
    To( L )=0;
    
    cnt=msk(im,:)+msk(ip,:)+msk(:,jm)+msk(:,jp);
    cnt( find(cnt==0) )=1;
    Tbr=(To(im,:)+To(ip,:)+To(:,jm)+To(:,jp))./cnt;
    
    To( L )=Tbr( L );
    
    To( find(To==0) )=NaN;
    
    T=To;
end
