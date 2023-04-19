function [out] = IntegrationMethods(a,b,f,n)
%%
% IntegrationMethods uses various methods to approximate an integral and
% reports the relative error.
%   [out] = IntegrationMethods(a,b,f,n) where
%       f is the anonymous function
%       a, b is the left and right endpoints of the interval a-b
%       out is the table indicating values of each method and the relative
%       error
%       n is the number of subdivisions (default n = 100)

%% initial values for testing
% a = -1; b = 1; n = 100;
% syms x; f = exp(1-x^3);

if nargin < 4 || isempty(n), n = 100; end
h = (b-a)/n;
xi = (a:h:b);
fprintf('Method \t Result \t Relative Error \t \n')

%% Actual integration
Actual = double(int(f,a,b));
RE_Actual = 0;
fprintf('Actual \t %.4f \t %.4f \n', Actual, RE_Actual);

%% Left endpoint UDF

%% Right endpoint UDF

%% Midpoint UDF

%% Trapezoid numerical integration (MATLAB built-in function 'trapz')
    function [int_trapz, err_trapz] = Int_Trapz(xi,f,Actual)
        f = @(x)f;
        y = f(xi);
        int_trapz = trapz(xi,y);
        err_trapz = abs((int_trapz-Actual)/Actual); % add equation for relative error
        fprintf('Trapezoid \t %.4f \t %.4f \n',int_trapz,err_trapz)
    end

%% Simpson_13 UDF (slightly modified from 'Numerical Methods for Engineers and Scientists Using MATLAB' 3e by Esfandiari) 
    function I_simp13 = Simpson_13(f,h,xi,n,Actual)
        f = @(x)f;
        I_simp13 = 0;
        for i = 1:2:n
            I_simp13 = I_simp13 + f(xi(i)) + 4*f(xi(i+1)) + f(xi(i+2));
        end
        I_simp13 = (h/3)*I_simp13;
        err_simp = abs((I_simp13-Actual)/Actual); % add equation for relative error
        fprintf('Simpson 1/3 \t %.4f \t %.4f \n',I_simp13,err_simp)
    end

%% Call all functions
Int_Trapz(xi,f);
Simpson_13(f,h,xi,n);

%% Summarized table output (to be deleted)
% fprintf('Method \t Result \t Relative Error \t \n')
% % LE
% % RE
% % MP
% % Trapezoid
% fprintf('Trapezoid \t %.f \t %.f \n',int_trapz,err_trapz)
% % Simpson 1/3

end

