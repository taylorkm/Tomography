% Create the square image, and get arrays ii and jj that have the relative
% coordinates.
nperdim = 128;
imageType = 'phantom';
if strcmp(imageType,'phantom')
    f = phantom(nperdim);
elseif strcmp(imageType,'square')
    xrng = linspace(-2, 2, nperdim);    
    yrng = linspace(-2, 2, nperdim);
    [x1,x2] = meshgrid(xrng, yrng);
    x = [x1(:), x2(:)];
    f = max(abs(x), [], 2) <= 1;% & abs(x(:,1))<.5;
end

figure
imagesc(reshape(f,nperdim,nperdim))
colormap gray
title('Original Image')


% coordinates of centers in x1x2-plane
[ii,jj] = meshgrid(1:nperdim, 1:nperdim);
ii = ii-ii(1, ceil(nperdim/2));
jj = jj-jj(ceil(nperdim/2), 1);
coords = [ii(:),jj(:)];
rads = sqrt(sum(coords.^2,2)); % radius of each coord from origin
maxrad = max(rads);

nTheta = 180;
thetas = pi*(0:nTheta-1)/nTheta;
nBins = ceil(1.45*nperdim);

bnds = - (ceil(maxrad)+1):2*(ceil(maxrad)+1)/nBins:(ceil(maxrad)+1);

R = zeros(nBins, nTheta);
for i = 1:nTheta
    th = thetas(i);
    % Build projector onto subspace orthogonal to direction of scans.
    u = [cos(th);sin(th)];     
    P = u*u'; % orthoprojector onto subspace spanned by u
    Q = [cos(th),-sin(th);sin(th),cos(th)];    
    PQ = P'*Q;       
    newCoords = coords*PQ;
    
    x = newCoords(:,1); % the only coordinates that matter
    [xsrtd,idx] = sort(x);
    
    % A crude quadrature
    j = 1;
    k = 1;
    while j <= nBins && k <= nperdim^2
        if ( xsrtd(k) < bnds(j+1) )
            R(j,i) = R(j,i) + f(idx(k));
            k = k+1;
        else            
            j = j+1;    
        end
    end
end
figure
imagesc(R)
colormap hot
title('Radon Transform')

% TODO: filter R before back projection

dtheta = abs(thetas(2)-thetas(1));

recon = zeros(nperdim, nperdim);

map2idx = @(x) round((x-bnds(1))/(bnds(2)-bnds(1))) + 1;

for t = 1:nTheta
    temp = ii*cos(thetas(t))+jj*sin(thetas(t));
    i0 = map2idx(temp);
    i0 = i0(:);
    recon = recon + reshape(R(i0,t),nperdim,nperdim);
end



figure
imagesc(recon)
colormap gray
title('Inverse Radon Transform')


