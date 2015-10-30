function [ cellImages, indOverlapList, indDelete ] = DC_ACD_removeOverlap( cellImages )
	indDelete = false(1,length(cellImages));
	indOverlapList = false(1,length(cellImages));
	nObjects = length(cellImages);
	for c1 = 1:nObjects-1
		for c2 = c1+1:nObjects
			if sum( (cellImages{c1}(:) > 0) & (cellImages{c2}(:) > 0) ) > 0
				indOverlapList(c1) = true;
				indOverlapList(c2) = true;
				
				%disp( [ 'Removing overlap between cell ' num2str(c1) ' and cell ' num2str(c2) ] );
				
				Io = (cellImages{c1} > 0) & (cellImages{c2} > 0);
				Is = ones(size(Io));
				Is( Io > 0 ) = 0;
				
				p1 = sum(Io(Io(:)>0)) / sum(cellImages{c1}(cellImages{c1}(:)>0));
				%disp(['  - Overlap is ' num2str(p1) ' % of cell ' num2str(c1) ]);
				
				p2 = sum(Io(Io(:)>0)) / sum(cellImages{c2}(cellImages{c2}(:)>0));
				%disp(['    Overlap is ' num2str(p2) ' % of cell ' num2str(c2) ]);
				
				if p1 < 0.1 && p2 < 0.1
					%disp([ '  - removing all of overlapping area..']);
					cellImages{c1}  = cellImages{c1} .* Is;
					cellImages{c2}  = cellImages{c2} .* Is;
				else
					b1 = sum(cellImages{c1}(:));
					b2 = sum(cellImages{c2}(:));
					%disp(['  - Brightness of cell ' num2str(c1) ' is ' num2str(b1) ]);
					%disp(['    Brightness of cell ' num2str(c2) ' is ' num2str(b2) ]);
					if ((b1 - b2) / b2) > 0.1
						%disp([' -> Giving overlapping area to brightest cell (' num2str(c1) ')']);
						%                         Is = imerode(Is, [1 1 1; 1 1 1; 1 1 1]);
						Is = imerode(Is, [0 1 0; 1 1 1; 0 1 0]);
						origSize = sum(sum(cellImages{c2}>0));
						imTemp  = cellImages{c2} .* Is;
						if sum(sum(imTemp>0)) < (0.5 * origSize)
							%disp(['    ... too little of cell ' num2str(c2) ' remains, scheduled for deletion...']);
							indDelete(c2) = true;
						else
							cellImages{c2} = imTemp;
						end
					elseif ((b2 - b1) / b1) > 0.1
						%disp([' -> Giving overlapping area to brightest cell (' num2str(c2) ')']);
						%                         Is = imerode(Is, [1 1 1; 1 1 1; 1 1 1]);
						Is = imerode(Is, [0 1 0; 1 1 1; 0 1 0]);
						origSize = sum(sum(cellImages{c1}>0));
						imTemp  = cellImages{c1} .* Is;
						if sum(sum(imTemp>0)) < (0.5 * origSize)
							%disp(['    ... too little of cell ' num2str(c1) ' remains, scheduled for deletion...']);
							indDelete(c1) = true;
						else
							cellImages{c1} = imTemp;
						end
					else
						%disp([' -> cells ' num2str(c1) ' and ' num2str(c1) ' can not be disentangled, please manually check ROI boundaries..']);
						indDelete(c1) = true;
						indDelete(c2) = true;
					end
				end
				
			end
		end
	end
end