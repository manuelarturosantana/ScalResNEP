function [val,inds] = inrectangle(val, xmin,xmax,ymin,ymax)
    % Filter values that are within the specified rectangle in the complex plane
    inds = real(val) > xmin & real(val) < xmax & imag(val) > ymin & imag(val) < ymax;
    val = val(inds);
end