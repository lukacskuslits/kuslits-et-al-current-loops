function [pcal_lat, pcal_lon] = calcCoordBounds(EstCoords, domx, domy)
    if EstCoords<0.5
       EstCoords=0.5;
    end
    pcal_lat = EstCoords(1:2:end,:);
    pcal_lon = EstCoords(2:2:end,:);
    pcal_lat(pcal_lat>domx) = domx;
    pcal_lon(pcal_lon>domy) = domy;
end