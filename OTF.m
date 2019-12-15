function H = OTF(ROWS,COLS)

% balls
% sigma1 = 0.000000000000009;
% sigma2 = 0.010000000001;

% keyboard
sigma1 = 0.001;
sigma2 = 0.00998;

ky = -64.0;
for r = 1:COLS
    kx = -36.0;
    for c = 1:ROWS
       OTF(r,c) = exp(-sigma1*kx*kx - sigma1*ky*ky) - exp(-sigma2*kx*kx - sigma2*ky*ky);
       kx = kx + .1 ;
    end
    ky = ky+.1;
end
max1 = max(OTF);
max2 = max(max1);
scale = 1.0/max2;
H = OTF.*scale;

