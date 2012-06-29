function average = ensemble_average(kT, tol)

if nargin < 2
    tol = 1e-14;
end

d = 1;
a = 2;
r0 = sqrt(2.0);
r = @(theta, phi) sqrt((cos(theta) - cos(phi)).^2 + (sin(theta) - sin(phi)).^2);
U = @(theta, phi) d * (1 - exp(-a * (r(theta, phi) - r0))).^2;

% The denominator function comes directly from Maple.
denominator = @(theta, phi) sqrt(cos(phi).^2-2.*cos(phi).^2.*cos(theta).^2+3+cos(theta).^2-2.*sin(theta).*sin(phi).*cos(theta).*cos(phi));

boltzmann_factor = @(theta, phi) exp(-U(theta, phi) ./ kT) ./ denominator(theta, phi);

partition_function = dblquad(boltzmann_factor, -pi, pi, -pi, pi, tol, @quadl);

integrate = @(h) dblquad(@(theta, phi) h(theta, phi) .* boltzmann_factor(theta, phi), ...
                             -pi, pi, -pi, pi, tol, @quadl) / partition_function;
    
average = integrate(r);
