%function runTripleChannelBuilder

strSesPath = 'D:\Data\Raw\imagingdata\20130612\';

strRedIm910 = 'reference\reference_t0_ch00.tif';
strGreIm910 = 'reference\reference_t0_ch01.tif';

for intSes=1:2
strRecPath = sprintf('xyt%02.0f%sxyt%02.0f_',intSes,filesep,intSes);
strRedIm810 = 't%05.0f_ch00.tif';
strGreIm810 = 't%05.0f_ch01.tif';

%get mean im
intUseIms = 1000;
matRed = zeros(512,512,intUseIms,'uint8');
matGreen = zeros(512,512,intUseIms,'uint8');
for intIm=1:intUseIms
	matRed(:,:,intIm) = imread([strSesPath strRecPath sprintf(strRedIm810,intIm)]);
	matGreen(:,:,intIm) = imread([strSesPath strRecPath sprintf(strGreIm810,intIm)]);
end
imRed810 = imnorm(mean(matRed,3)/255);
imRed910 = im2double(imread([strSesPath strRedIm910]));
imGreen910 = im2double(imread([strSesPath strGreIm910]));
imGreen810 = imnorm(mean(matGreen,3)/255);

%register
output = dftregistration(fft2(imRed810),fft2(imRed910),20);
dblYCorrect = output(3);
dblXCorrect = output(4);

imRed910 = imnorm(imenhance(circshift(imRed910,[round(dblYCorrect) round(dblXCorrect)])));
imGreen910 = imnorm(imenhance(circshift(imGreen910,[round(dblYCorrect) round(dblXCorrect)])));
imRed810 = imnorm(imenhance(imRed810));
imGreenChan = imnorm(imenhance(imGreen810));

%rectify & normalize
imRedChan = imRed810;%./(imRed910+imRed810);
%imRedChan(imRedChan<0) = 0;
imRedChan = imnorm(imRedChan);
imBlueChan = imRed910;%./(imRed910+imRed810);
%imBlueChan(imBlueChan<0) = 0;
imBlueChan = imnorm(imBlueChan);
%imBlueChan = zeros(512,512);

%assign to RGB
imRGB = zeros(512,512,3);
imRGB(:,:,1) = imRedChan;
imRGB(:,:,2) = imGreenChan;
imRGB(:,:,3) = imBlueChan;

figure

imshow(imRGB);
draw
end