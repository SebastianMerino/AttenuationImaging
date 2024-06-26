Bmode1 = getBmode(rf(:,:,1));
Bmode2 = getBmode(rf(:,:,2));
Bmode3 = getBmode(rf(:,:,3));
dynRange = [-50 0];

%%
Bmode = getBmode(mean(rf,3));


figure,tiledlayout(1,4)
nexttile,
imagesc(x,z,Bmode1, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,Bmode2, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,Bmode3, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,Bmode, dynRange)
colormap gray
axis image


%%
function bm = getBmode(rf)
    bm = db(hilbert(rf));
    bm = bm - max(bm(:));
end

