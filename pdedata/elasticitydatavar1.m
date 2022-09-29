function pde = elasticitydatavar1(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

para = struct('lambda',lambda,'mu',mu);
[u1,u2,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(para);

% rhs
    function val = f(p)  
        x = p(:,1); y = p(:,2);
        val = [f1(x,y), f2(x,y)];
    end
% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = u1(x,y);
        val(:,2) = u2(x,y);
    end
% Dirichlet
    function val = g_D(p)
        val = uexact(p);
    end
% derivative of the exact solution
    function val =  Du(p)  % u = [u1,u2]
        x = p(:,1); y = p(:,2);
        val(:,1) = u1x(x,y);  val(:,2) = u1y(x,y); 
        val(:,3) = u2x(x,y);  val(:,4) = u2y(x,y);
    end

pde = struct('lambda',lambda, 'mu',mu, 'f', @f, 'exactu',@uexact,'g_D',@g_D,'Du',@Du);
end

function [u1,u2,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(para)    
    lambda = para.lambda; mu = para.mu; 
    syms x y;    
    
    % exact solution
    u1 = (-1+cos(2*pi*x))*sin(2*pi*y) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);
    u2 = -(-1+cos(2*pi*y))*sin(2*pi*x) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);

    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);
    u2x = diff(u2,x);      u2y = diff(u2,y);

    % Laplace u
    Lap1 = diff(u1,x,2)+diff(u1,y,2);
    Lap2 = diff(u2,x,2)+diff(u2,y,2);

    % grad (div u)
    divu = u1x+u2y;
    graddiv1 = diff(divu,x);
    graddiv2 = diff(divu,y);

    % Linear elasticity
    Ls1 = mu*Lap1 + (lambda+mu)*graddiv1;
    Ls2 = mu*Lap2 + (lambda+mu)*graddiv2;

    % f
    f1 = - Ls1; %f1 = simplify(f1);
    f2 = - Ls2; %f2 = simplify(f2);
    
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y});    u2 = matlabFunction(u2,'Vars',{x,y});
    f1 = matlabFunction(f1,'Vars',{x,y});    f2 = matlabFunction(f2,'Vars',{x,y}); 
    u1x = matlabFunction(u1x,'Vars',{x,y});  u1y = matlabFunction(u1y,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});  u2y = matlabFunction(u2y,'Vars',{x,y});     
end


