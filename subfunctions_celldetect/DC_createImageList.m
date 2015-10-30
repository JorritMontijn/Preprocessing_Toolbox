function [cellIm,cellImName] = DC_createImageList(imProc)
	%DC_createImageList Creates image list from 2 image stacks
	%   [cellIm,cellImName] = DC_createImageList(imRGB,imHQ)
	
	%list; 1=normal overlay; 2=HQ overlay; 3+ =single channels (normal
	%first)
	
	%normalize image
	imRGB = zeros(size(imProc));
	for intCh=1:size(imProc,3)
		imRGB(:,:,intCh) = imadjust(imProc(:,:,intCh));
	end
	imProc(isnan(imProc)) = 0;
	
	%create enhanced image
	imHQ = imenhance(imProc);
	
	%im 1
	cellIm{1} = imRGB;
	cellImName{1} = 'Overlay normalized images';
	
	%im 2
	cellIm{2} = imHQ;
	cellImName{2} = 'Overlay contrast-enhanced images';
	
	%loop through rest and add to cells if non-zero
	im3Ch = zeros(size(imRGB));
	intCh = 2;
	for intSingleChannel=1:6
		%get channel
		intIm = mod(intSingleChannel,3);
		if intIm == 0, intIm = 3;end
		
		%decide which stack it's from
		if intSingleChannel < 4
			imStack = imRGB;
			strType = 'normalized';
		else
			imStack = imHQ;
			strType = 'contrast-enhanced';
		end
		
		%define image from channel
		imThis = imStack(:,:,intIm);
		
		%get min/max vals
		imMin = min(imThis(:));
		imMax = max(imThis(:));

		%check if non-zero
		if imMin ~= imMax
			%increment channel
			intCh = intCh + 1;

			%create overlay
			imThis3 = im3Ch;
			imThis3(:,:,intIm) = imThis;
			
			%add to cells
			cellIm{intCh} = imThis3; %#ok<AGROW>
			cellImName{intCh} = sprintf('Channel %d %s',intIm,strType); %#ok<AGROW>
		end
	end
end

