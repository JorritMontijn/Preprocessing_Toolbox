function imI = DC_ACD_thresholdContour( imI, dblPercentage )
    % get threshold value
    ImV = imI( imI(:) > 0 );
    ImV = sort( ImV(:) ) ;
    dblThresh = ImV( round(length(ImV) * dblPercentage) ) ;
    
    % threshold image
    imI( imI < dblThresh ) = 0 ;
end