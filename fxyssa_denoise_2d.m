function [ D1 ] = fxyssa_denoise_2d(D,flow,fhigh,dt,N)
%  FXY_SSA: F-XY domain singular spectrum analysis (SSA)
%  Sacchi, M. D. (2009). FX singular spectrum analysis. CSPG CSEG CWLS Convention.
	
if nargin==2
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
end

[nt,nx,ny]=size(D);

nf=2^nextpow2(nt);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,nx,ny);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1
    ilow=1;
end

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1
    ihigh=floor(nf/2)+1;
end

lx=floor(nx/2)+1;
ly=floor(ny/2)+1;

% main loop
for k=ilow:ihigh
    
    if(ny==1)
        S_obs=DATA_FX(k,:,:).';
    else
        S_obs=squeeze(DATA_FX(k,:,:));   
    end
    Sn_1=S_obs;
        
    M=P_H(Sn_1,lx,ly);  % form hankel
    M=P_R(M,N);         % compute svd
    Sn=P_A(M,nx,ny,lx,ly);  % average
    
    for j=1:ny
        DATA_FX0(k,:,j) = DATA_FX0(k,:,j)+reshape(Sn(:,j),1,nx);
    end
    
end

% symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:,:) = conj(DATA_FX0(nf-k+2,:,:));
end

% Back to TX
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:nt,:,:);

return

function [dout]=P_H(din,lx,ly)
% forming block Hankel matrix
[nx,ny]=size(din);
lxx=nx-lx+1;

for j=1:ny
    r=hankel(din(1:lx,j),[din(lx:nx,j)]);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r;
        end
    else
        for id=1:(ny-j+1)
            dout((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r;
        end
    end
end
return

function [dout]=P_R(din,N)
[U,D,V]=svd(din);
dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');
return

function [dout]=P_A(din,nx,ny,lx,ly)
% Averaging the block Hankel matrix to output the result
lxx=nx-lx+1;
dout=zeros(nx,ny);

for j=1:ny
    if j<ly
        for id=1:j
            dout(:,j) =dout(:,j)+ ave_antid(din(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
        end
    else
        for id=1:(ny-j+1)
            dout(:,j) =dout(:,j)+ ave_antid(din((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(ny-j+1);
        end
    end
end
return


function [dout] =ave_antid(din)
% averaging along antidiagonals
[n1,n2]=size(din);
nout=n1+n2-1;
dout=zeros(nout,1);
for i=1:nout
    if i<n1
        for id=1:i
            dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i;
        end
    else
        for id=1:nout+1-i
            dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i);
        end
    end
end
return
