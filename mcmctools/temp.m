[xp,yp] = deal(0:.02:1);
[x,y] = meshgrid(xp,yp);
f = exp(x+y) + sin((x-2*y)*3);
fn = f + randn(size(f))*0.5; % Adding white noise
fs = smoothn(fn); % DCT-based smoothing
subplot(121), surf(xp,yp,fn), zlim([0 8]), axis square
subplot(122), surf(xp,yp,fs), zlim([0 8]), axis square