function h = atmospalt( P )
%  ATMOSPALT Calculate altitude based on ambient pressure.

if ~isnumeric( P )
    % Altitude should be a numeric array.  Otherwise error.
    error('aero:atmospalt:notnumeric', ...
        'The pressure altitude input must be a numeric value.');
end

%
% From https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf
h = (1-(P./1013.25).^0.190284 ).*145366.45.*0.3048;

%From https://en.wikipedia.org/wiki/Barometric_formula
% h = -8400.*log(P/1013.25);
% p = p0e-(h/h0)

end
