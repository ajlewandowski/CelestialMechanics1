% This function convert the r and v state vectors to the classical orbital
% elements. 
%
% output is semi-major axis, eccentricity, inclination, arg of periapsis, 
% ascending node, mean anomaly in units of degrees
%
% Alec Lewandowski
% 10/6/25

function [oe] = eci2oeNU(r, v, mu)

 h = cross(r,v);
    hnorm = norm(h);
    
    % Calculate the node vector
    n = cross([0 0 1], h);
    nmag = norm(n);
    
    % Calculate the eccentricity vector
    evec = (1/mu)*((norm(v)^2 - mu/norm(r)) * r - dot(r,v) * v);
    e = norm(evec);
    
    % Calculate the inclination (i)
    inc = acos(h(3)/hnorm);
    
    % Calculate RAAN (Ω)
    if nmag ~= 0
        RAAN = atan2(n(2), n(1)); % Right Ascending Node
    else
        RAAN = 0;
    end
    
    % Calculate Argument of Perigee (ω)
    if e > 1e-12 && nmag ~= 0
        argp = atan2(dot(cross(n, evec), h)/hnorm, dot(n, evec)); % Argument of Perigee
    else
        argp = 0;
    end
    
    % Calculate True Anomaly (ν)
    if e > 1e-12
        nu = atan2(dot(cross(evec, r), h)/hnorm, dot(evec, r)); % True anomaly
    else
        nu = atan2(dot(cross(n, r), h)/hnorm, dot(n, r)); % For circular orbit
    end
    
    % Semi-major axis (a) using the vis-viva equation
    sma = 1 / (2/norm(r) - norm(v)^2/mu);
    
    % Ensure the true anomaly is between 0 and 2π
    nu = mod(nu, 2*pi);
    
    % Return orbital elements: [semi-major axis, eccentricity, inclination, argument of periapsis, RAAN, true anomaly]
    oe = [sma, e, inc, argp, RAAN, nu];
end