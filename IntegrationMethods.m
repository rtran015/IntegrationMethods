function [out] = IntegrationMethods(a,b,f,n)
%%
% IntegrationMethods uses various methods to approximate an integral and
% reports the relative error.
%   [out] = IntegrationMethods(a,b,f,n) where
%       f is the symbolic function
%       a, b is the left and right endpoints of the interval a-b
%       out is the table indicating values of each method and the relative
%       error
%       n is the number of subdivisions (default n = 100)

%% initial values for testing
% a = -1; b = 1; n = 100; h = (b-a)/n; xi = (a:h:b); syms x; f = exp(1-x^3); ff = matlabFunction(f);

%% Initialize variables
if nargin < 4 || isempty(n), n = 100; end
h = (b-a)/n;
xi = (a:h:b);
ff = matlabFunction(f);

% Actual integration
Actual = double(int(f,a,b));
RE_Actual = 0;

%% Output
fprintf('Integration of f = %s, interval [%i,%i], n = %i\n',f,a,b,n);
fprintf('Method \t\t\t\t Result \t\t\t Relative Error \t \n')
fprintf('Actual \t\t\t\t %.4f \t\t\t %.4f%% \n', Actual, RE_Actual);

% Call all functions
Int_LE(xi,ff,Actual);
Int_RE(xi,ff,Actual);
Int_MP(xi,ff,Actual);
Int_Trapz(xi,ff,Actual);
Simpson_13(ff,h,xi,n,Actual);

%% Left endpoint UDF
    function [int_LE, err_LE] = Int_LE(xi,ff,Actual)
        int_LE = 0;
        for i = 1:1:n
            int_LE = int_LE + h*ff(xi(i));
        end
        err_LE = abs((int_LE-Actual)/Actual)*100;
        fprintf('Left Endpoint \t\t %.4f \t\t\t %.4f%% \n',int_LE,err_LE)
    end

%% Right endpoint UDF
    function [int_RE, err_RE] = Int_RE(xi,ff,Actual)
        int_RE = 0;
        for i = 2:1:n+1
            int_RE = int_RE + h*ff(xi(i));
        end
        err_RE = abs((int_RE-Actual)/Actual)*100;
        fprintf('Right Endpoint \t\t %.4f \t\t\t %.4f%% \n',int_RE,err_RE)
    end

%% Midpoint UDF
    function [int_MP, err_MP] = Int_MP(xi,ff,Actual)
        int_MP = 0;
        for i = 1:1:n
            int_MP = int_MP + h*ff((xi(i)+xi(i+1))/2);
        end
        err_MP = abs((int_MP-Actual)/Actual)*100;
        fprintf('Midpoint \t\t\t %.4f \t\t\t %.4f%% \n',int_MP,err_MP)
    end

%% Trapezoid numerical integration (MATLAB built-in function 'trapz')
    function [int_trapz, err_trapz] = Int_Trapz(xi,ff,Actual)
        % ff = matlabFunction(f);
        y = ff(xi);
        int_trapz = trapz(xi,y);
        err_trapz = abs((int_trapz-Actual)/Actual)*100;
        fprintf('Trapezoid \t\t\t %.4f \t\t\t %.4f%% \n',int_trapz,err_trapz)
    end

%% Simpson_13 UDF (slightly modified from 'Numerical Methods for Engineers and Scientists Using MATLAB' 3e by Esfandiari) 
    function I_simp13 = Simpson_13(ff,h,xi,n,Actual)
        % ff = matlabFunction(f);
        I_simp13 = 0;
        for i = 1:2:n
            I_simp13 = I_simp13 + ff(xi(i)) + 4*ff(xi(i+1)) + ff(xi(i+2));
        end
        I_simp13 = (h/3)*I_simp13;
        err_simp = abs((I_simp13-Actual)/Actual)*100;
        fprintf('Simpson 1/3 \t\t %.4f \t\t\t %.4f%% \n',I_simp13,err_simp)
    end
end

