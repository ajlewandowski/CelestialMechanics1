% This function runs solution for HW4 Problem 3
%
% Alec Lewandowski
% 10/22/25
clear variables
close all
clc

format long
% Defining conditions given
phi0 = deg2rad(62);
n = 10;

%% Problem 3.i
x = sin(phi0);

syms t
P = sym(zeros(1, n+1));  
P(1) = sym(1);            
P(2) = t;                  

for i = 1:n-1
    P(i+2) = simplify(((2*i+1)*t*P(i+1) - i*P(i)) / (i+1));
end

Pi_values = double(subs(P, t, x)); 

disp('Pi(sin(phi)) for i=0 to n:');
disp(Pi_values');

P_ij = sym(zeros(n));  
for i = 1:n
    for j = 1:i
        var1 = 1/(2^i*factorial(i));
        var2 = (1-t^2)^(j/2);
        var3 = diff((t^2 - 1)^i,t,i+j);
        P_ij(i,j) = var1*var2*var3;
    end
end

Pij_final = double(subs(P_ij, t, x)); 

disp('Matrix P_ij:');
disp(Pij_final);

Pi_norm = norm(Pi_values, 2);     
Pij_norm = norm(Pij_final, 2);        

disp(['2-norm of Pi vector: ', num2str(Pi_norm)]);
disp(['2-norm of P_ij matrix: ', num2str(Pij_norm)]);

%% Problem 3.ii
[J, C, S, Rstar, GMstar] = getMoonData();

% First going to compute the gravity potentials
lat = deg2rad(-89.9999:3:89.9999);
lon = deg2rad(-179.999:3:179.999);

V = zeros(length(lat), length(lon));

r = Rstar; % because we are evaluating surface

for i = 1:length(lat)
    for j = 1:length(lon)
        phi = (lat(i));
        lambda = (lon(j));         
        potential = 0;                 
        U1 = 0;
        % Loop over degrees l and orders m
        for l = 1:20 
            Pl = legendre(l, sin(phi), 'norm');

            temp = ((Rstar/r)^l) * Pl(l)*J(l);

            for m = 1:l
                
                P_ij = sym(zeros(n));  
                for b = 1:n
                    for c = 1:i
                        var1 = 1/(2^b*factorial(b));
                        var2 = (1-t^2)^(c/2);
                        var3 = diff((t^2 - 1)^b,t,b+c);
                        P_ij(b,c) = var1*var2*var3;
                    end
                end
                
                Pij_final = double(subs(P_ij, t, x));
                
                if m == 1
                    Clm = J(l);      
                    Slm = 0;         
                else
                    Clm = C(l, m);
                    Slm = S(l, m);
                end
                term = Plm(l+m) * (Clm * cos(m * lambda) + Slm * sin(m * lambda));
                potential = potential + (Rstar / r)^l * term;
            end

            U1 = U1 + temp;

        end

        V(i, j) = (GMstar/r)*(-U1 + potential);
    end
end

figure;
surf(rad2deg(lon), rad2deg(lat), V);  % Create a 3D surface plot
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
zlabel('Lunar Potential (m^2/s^2)');
title('Lunar Non-Spherical Gravitational Potential');
view(0,90)
colorbar;