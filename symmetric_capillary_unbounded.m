% function [] = symmetric_capillary_unbounded()
function [sigma, Tvec, sig_min, sig_max, bvec ] = symmetric_capillary_unbounded()
% Ray Treinen, June 2022
%
% Compute various quantities related to the capillary surfaces known as
% symmetric unbounded liquid bridges.  These surfaces have generating
% curves that contain vertial points at a point (sigma,T(sigma)) and are
% asypmtotic to the horizontal axis.
%
% The output can be set to [sigma, Tvec, sig_min, sig_max, bvec ]
% to return the solution to the values of sigma and T(sigma)
% The default is set to merely plot graphs of interest.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.


%% physical parameters

sig_min = 0.085;
sig_max = 2;

kappa = 1;
a = 1;
b = max(14,a + 4);
psia = -pi/2;
psib = 0;

%tic
%% Computational parameters
% computational domain
X = [-1,1];

% maximum loop counts and the tolerances for each loop
max_iter_newton = 100000;
max_iter_bvp = 10000;
max_abs = 20;
tol_newton = 1e-14;
tol_bvp = 1e-13;
tol_abs = 1e-11;
options = odeset('AbsTol',1e-11,'RelTol',1e-11);

% initialize the number of points.
k = 7;
n = 2*k + 1;

%% Determining the value of T(sigma) at sigma = a
bigN = 100;
sigma = chebpts(bigN,[sig_min;sig_max]);
Tvec = zeros(size(sigma));
bvec = Tvec;

for ii = 1:bigN
    
    a = sigma(ii);
    b = max(14,a + 4);
    
    k = 7;
    n = 2*k + 1;
    s = chebpts(n);

    R0 = @(s) (1 + s)*b/2 + (1 - s)*a/2;
    U0 = @(s) exp(-R0(s) + a);
    Psi0 = @(s) atan(-exp(-R0(s) + a));
    ell0 = b - a;
    v = [R0(s); U0(s); Psi0(s); ell0];
    
    res_abs = 1;
    iter_abs = 0;
    T = sqrt(2);
    while res_abs > tol_abs
        [v, n] = cheb_engine(v, n);
        R = v(1:n);
        U = v(n+1:2*n);
        %     figure(14)
        %     plot(chebfun(R),chebfun(U),'.-k')
        %     hold on
        %     drawnow
        %     Psi = v(2*n+1:end-1);
        %     ell = v(end);
        res_abs = abs(T - U(1));
        if res_abs > tol_abs
            R(R>1) = R(R>1)*(b + 2)/b;
            b = b + 2;
            v(1:n) = R;
        end
        T = U(1);
        iter_abs = iter_abs + 1;
        if iter_abs > max_abs
            disp('Maximum number of iterations reached in finding T')
            break
        end
    end
%     figure(1)
%     plot(chebfun(R),chebfun(U),'.-k')
%     axis([0 b 0 sqrt(2)])
%     hold on
%     
%     
    Tvec(ii) = T;
    bvec(ii) = b;
end


% figure(2)
% plot(chebfun(sigma),chebfun(Tvec),'.-k')
% axis([0 sig_max 0 sqrt(2)])
% title('T(\sigma)')
% hold on
% sg = chebpts(100,[sig_min/100; 2*sig_min]);
% Tsg = -sg.*log(sg);
% plot(chebfun(sg),chebfun(Tsg),'k')
% 
% siginterp = chebpts(100,[sig_min/100; sig_min]);
% Tinterp = interp1([sg(1:10);sigma],[Tsg(1:10);Tvec],siginterp,'spline');
% plot(chebfun(siginterp),chebfun(Tinterp),'--k')

DT = diffmat(bigN, 1, [sig_min;sig_max]);
Tprime = DT*Tvec;

% figure(21)
% plot(chebfun(sigma),chebfun(Tprime),'.-k')
% title('T^\prime(\sigma)')
% hold on
% % DTsg = -log(sg) - 1;
% % plot(chebfun(sg),chebfun(DTsg),'k')
% % 
% % DTinterp = interp1([sg(1:10);sigma],[DTsg(1:10);Tprime],siginterp,'pchip');
% % plot(chebfun(siginterp),chebfun(DTinterp),'--k')
% axis equal

return

%% Solving the ODE system for /dot r

rdotmin = zeros(size(Tvec));
rdotend = zeros(size(Tvec));

for ii = 1:length(Tvec)
    sig = sigma(ii);
    Tsig = Tvec(ii);
    DTsig = Tprime(ii);
    
    ic = [sig;Tsig;1;DTsig];
    [phi,vvv] = ode45(@bigODE,[pi/2 0],ic,options);

    rode = vvv(:,1);
    uode = vvv(:,2);
    rdot = vvv(:,3);
    udot = vvv(:,4);
    rdotmin(ii) = min(rdot);
    rdotend(ii) = rdot(end);
    
%     figure(2)
%     plot(rode,uode,'k')
%     axis equal
    
    
    figure(6)
    plot(phi,rdot,'k')
    axis equal
    hold on
    
    % pause
end
title('Foliating a portion of the \phi rdot-plane')

figure(7)
plot(sigma,rdotmin,':k')
hold on
plot(sigma,rdotend,'k')
title('Select values of rdot as a function of \sigma')

% %% pass through the interpolated data
% Uncomment the code after plotting figures 2 and 21.
% 
% rdotmin = zeros(size(Tinterp));
% rdotend = zeros(size(Tinterp));
% 
% for ii = 1:length(Tinterp)
%     sig = siginterp(ii);
%     Tsig = Tinterp(ii);
%     DTsig = DTinterp(ii);
% %     DTsig(sig < 1) = -log(sig(sig < 1)) - 1;
%     ic = [sig;Tsig;1;DTsig];
%     [phi,vvv] = ode45(@bigODE,[pi/2 0],ic,options);
% 
%     rode = vvv(:,1);
%     uode = vvv(:,2);
%     rdot = vvv(:,3);
%     udot = vvv(:,4);
%     rdotmin(ii) = min(rdot);
%     rdotend(ii) = rdot(end);
%     
% %     figure(2)
% %     plot(rode,uode,'k')
% %     axis equal
%     
%     
%     figure(8)
%     plot(phi,rdot,'k')
%     axis equal
%     hold on
%     
%     % pause
% end
% 
% figure(9)
% plot(siginterp,rdotmin,':k')
% hold on
% plot(siginterp,rdotend,'k')

%% the main computational engine of the file
    function [v, n] = cheb_engine(v, n)
        
        % intialize the residual
        res_bvp = 1;
        
        while res_bvp > tol_bvp
            
            % initialize the differential operator components
            %
            % D0 and D1 are spectral operators corresponding to a
            % downsampling matrix and a differentiation matrix,
            % respectively
            D0 = diffmat([n-1 n],0,X);
            D1 = diffmat([n-1 n],1,X);
            Z0 = sparse(n-1,n);
            D01 = spalloc(n-1, 3*n + 1, n*(n-1));
            D02 = D01;
            D03 = D01;
            D11 = D01;
            D12 = D01;
            D13 = D01;
            D01(1:n-1, 1:n) = D0;
            D11(1:n-1, 1:n) = D1;
            D02(1:n-1, n+1:2*n) = D0;
            D12(1:n-1, n+1:2*n) = D1;
            D03(1:n-1, 2*n+1:3*n) = D0;
            D13(1:n-1, 2*n+1:3*n) = D1;
            
            % Evaluate the computational vector to check the boundary
            % conditions
            dT1n1 = sparse(1,3*n+1);
            dT1p1 = dT1n1;
            dT3n1 = dT1n1;
            dT3p1 = dT1n1;
            dT1n1(1) = 1;
            dT1p1(n) = 1;
            dT3n1(2*n+1) = 1;
            dT3p1(end -1) = 1;
            
            % building the nonlinear operator N
            % and the linear operator L
            
            if sqrt(kappa)*a <= 1
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    (D01*v).*(D13*v) + v(end).*sin(D03*v) - kappa*v(end).*(D02*v).*(D01*v)
                    dT1n1*v - a
                    dT1p1*v - b
                    dT3n1*v - psia
                    dT3p1*v - psib ];
                
                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(D13*v - kappa*v(end).*(D02*v),0,n-1,n-1)*D0, spdiags(-kappa*v(end)*(D01*v),0,n-1,n-1)*D0, spdiags(D01*v,0,n-1,n-1)*D1 + spdiags(v(end)*cos(D03*v),0,n-1,n-1)*D0, sin(D03*v) - kappa*(D02*v).*(D01*v)
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            else
                
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    D13*v + v(end).*sin(D03*v)./(D01*v) - kappa*v(end).*D02*v
                    dT1n1*v - a
                    dT1p1*v - b
                    dT3n1*v - psia
                    dT3p1*v - psib ];
                
                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(-v(end)*sin(D03*v)./((D01*v).^2),0,n-1,n-1)*D0, -kappa*v(end)*D0, D1 + (spdiags(v(end)*cos(D03*v),0,n-1,n-1))*D0, sin(D03*v)./(D01*v) - kappa*D02*v
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            end
            
            % initialize a counter and
            % the residual for the Newton's method loop
            kk = 1;
            res_newton = 1;
            
            %% Newton's method loop
            while res_newton > tol_newton
                
                lastwarn('', '');
                % the main steps
                dv = -L(v)\N(v);
                
                % warning syntax to catch badly scaled matrices.  This
                % happens when b is too small.
                [warnMsg, warnID] = lastwarn();
                if(isempty(warnID))
                else
                    warning('Radii and inclination angles lead to a muilt-scale problem.')
                    return
                    % plot the current configuration if there is a problem.
                    %                    R10 = v(1:n)
                    %                     U10 = v(n+1:2*n)
                    %                     Psi10 = v(2*n+1:end-1)
                    %                     ell10 = v(end)
                    %                     figure(12)
                    %                     plot(chebfun(R10),chebfun(U10),'.-k')
                    %
                    %                     pause
                    %
                    %                     temp_psi = interp1(X,[psia;psib],chebpts(length(2*n+1:3*n)));
                    %                     v(2*n+1:3*n) = temp_psi;
                end
                
                % Newton's step
                v = v + dv;
                
                % the barrier if the solution strays too far from the
                % expected solution.  This forces the inclination angle of
                % the solution to remain within reasonable bounds
                temp_psi = v(2*n+1:3*n);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n+1:3*n) = temp_psi;
                
                temp_R = v(1:n);
                temp_R(temp_R <= 0) = a/2;
                v(1:n) = temp_R;
                
                if v(end) <= (b - a)/2
                    v(end) = b-a;
                elseif v(end) > pi*(b + a)
                    v(end) = ell0;
                end
                
                % computing the residual of the Newton's method
                res_newton = norm(dv,'fro')/norm(v,'fro');
                
                % updating the counter and if it grows above the prescribed
                % mazimum, stop this loop and report to the outer loop.
                kk = kk+1;
                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
                    break
                end
                
            end
            
            %% residual and other loop conditions
            
            % the relative norm of the residual
            res_bvp = norm(N(v),'fro')/norm(v,'fro');
            
            % adaptively add more Chebyshev points and resample the state
            % of the problem if the tolerance is not met.
            % Additionally, if there is too much numerical oscillation, add
            % more Chebyshev points and resample the state
            % of the problem and reset the residual.
            S2 = diffmat([2*n n],0,X);
            if res_bvp > tol_bvp
                nn = n;
                n = n + 4;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;
                
            elseif length(find(diff(sign(diff((S2*v(2*n+1:3*n))))))) >= 2
                nn = n;
                %                 n = 2*(n - 1) - 1;
                n = n + 4;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;
                
                res_bvp = 1;
            else
                break
            end
            
            % if the function exceeds the maximum number of iterations,
            % break with an error statement.
            if n > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                break
            end
            
            
        end
    end

%% ODE function

    function out = bigODE(phi,vv)
        out = [-vv(1).*cos(phi)./(vv(1).*vv(2) + sin(phi)); ...
            -vv(1).*sin(phi)./(vv(1).*vv(2) + sin(phi)); ...
            cos(phi).*(vv(4).*vv(1).^2 - vv(3).*sin(phi))./(vv(1).*vv(2) + sin(phi)).^2; ...
            sin(phi).*(vv(4).*vv(1).^2 - vv(3).*sin(phi))./(vv(1).*vv(2) + sin(phi)).^2 ];
    end

end