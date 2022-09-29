function w = Base2d_P3(wStr,node,elem,quadOrder)
%      3
%      |  \
%      5    7
%      |     \
%      8  10  4
%      |        \
%      1- 6 - 9 - 2
% The order of the d.o.f.s is shown above.

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder);  
ng = length(weight);  
NT = size(elem,1);

dot = strfind(wStr,'.');
wStr = wStr(dot:end);

% gradbasis
if ~strcmpi(wStr,'.val')
    Dlambda = gradbasis(node,elem);
    Dlambda1 = Dlambda(1:NT,:,1);
    Dlambda2 = Dlambda(1:NT,:,2);
    Dlambda3 = Dlambda(1:NT,:,3);
    Dlambdax = [Dlambda1(:,1), Dlambda2(:,1), Dlambda3(:,1)];
    Dlambday = [Dlambda1(:,2), Dlambda2(:,2), Dlambda3(:,2)];
end

%% u.val
if strcmpi(wStr,'.val')
    w1 = 1/2*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
    w2 = 1/2*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
    w3 = 1/2*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
    w4 = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,2)-1);
    w5 = 9/2*lambda(:,3).*lambda(:,1).*(3*lambda(:,3)-1);
    w6 = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
    w7 = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,3)-1);
    w8 = 9/2*lambda(:,3).*lambda(:,1).*(3*lambda(:,1)-1);
    w9 = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);
    w10 = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);
    w1 = repmat(w1',NT,1); w2 = repmat(w2',NT,1); w3 = repmat(w3',NT,1);
    w4 = repmat(w4',NT,1); w5 = repmat(w5',NT,1); w6 = repmat(w6',NT,1);
    w7 = repmat(w7',NT,1); w8 = repmat(w8',NT,1); w9 = repmat(w9',NT,1);
    w10 = repmat(w10',NT,1);
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dx
if strcmpi(wStr,'.dx')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = 1/2*(27*lambda(p,1)^2-18*lambda(p,1)+2)*Dlambdax(:,1);
        w2(:,p) = 1/2*(27*lambda(p,2)^2-18*lambda(p,2)+2)*Dlambdax(:,2);
        w3(:,p) = 1/2*(27*lambda(p,3)^2-18*lambda(p,3)+2)*Dlambdax(:,3);
        w4(:,p) = 9/2*Dlambdax(:,3).*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*Dlambdax(:,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambdax(:,2));
        w5(:,p) = 9/2*Dlambdax(:,1).*lambda(p,3).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*Dlambdax(:,3).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambdax(:,3));
        w6(:,p) = 9/2*Dlambdax(:,1).*lambda(p,2).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambdax(:,2).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambdax(:,1));
        w7(:,p) = 9/2*Dlambdax(:,3).*lambda(p,2).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*Dlambdax(:,2).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambdax(:,3));
        w8(:,p) = 9/2*Dlambdax(:,1).*lambda(p,3).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambdax(:,3).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambdax(:,1));
        w9(:,p) = 9/2*Dlambdax(:,1).*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*Dlambdax(:,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambdax(:,2));
        w10(:,p) = 27*Dlambdax(:,1).*lambda(p,2).*lambda(p,3) ...
            + 27*lambda(p,1).*Dlambdax(:,2).*lambda(p,3) ...
            + 27*lambda(p,1).*lambda(p,2).*Dlambdax(:,3);
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dy
if strcmpi(wStr,'.dy')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = 1/2*(27*lambda(p,1)^2-18*lambda(p,1)+2)*Dlambday(:,1);
        w2(:,p) = 1/2*(27*lambda(p,2)^2-18*lambda(p,2)+2)*Dlambday(:,2);
        w3(:,p) = 1/2*(27*lambda(p,3)^2-18*lambda(p,3)+2)*Dlambday(:,3);
        w4(:,p) = 9/2*Dlambday(:,3).*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*Dlambday(:,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambday(:,2));
        w5(:,p) = 9/2*Dlambday(:,1).*lambda(p,3).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*Dlambday(:,3).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambday(:,3));
        w6(:,p) = 9/2*Dlambday(:,1).*lambda(p,2).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambday(:,2).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambday(:,1));
        w7(:,p) = 9/2*Dlambday(:,3).*lambda(p,2).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*Dlambday(:,2).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambday(:,3));
        w8(:,p) = 9/2*Dlambday(:,1).*lambda(p,3).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambday(:,3).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambday(:,1));
        w9(:,p) = 9/2*Dlambday(:,1).*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*Dlambday(:,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambday(:,2));
        w10(:,p) = 27*Dlambday(:,1).*lambda(p,2).*lambda(p,3) ...
            + 27*lambda(p,1).*Dlambday(:,2).*lambda(p,3) ...
            + 27*lambda(p,1).*lambda(p,2).*Dlambday(:,3);
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.grad
if strcmpi(wStr,'.grad')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,2*p-1:2*p) = 1/2*(3*Dlambda1).*(3*lambda(p,1)-2).*lambda(p,1) ...
            + 1/2*(3*lambda(p,1)-1).*(3*Dlambda1).*lambda(p,1) ...
            + 1/2*(3*lambda(p,1)-1).*(3*lambda(p,1)-2).*Dlambda1;
        w2(:,2*p-1:2*p) = 1/2*(3*Dlambda2).*(3*lambda(p,2)-2).*lambda(p,2) ...
            + 1/2*(3*lambda(p,2)-1).*(3*Dlambda2).*lambda(p,2) ...
            + 1/2*(3*lambda(p,2)-1).*(3*lambda(p,2)-2).*Dlambda2;
        w3(:,2*p-1:2*p) = 1/2*(3*Dlambda3).*(3*lambda(p,3)-2).*lambda(p,3) ...
            + 1/2*(3*lambda(p,3)-1).*(3*Dlambda3).*lambda(p,3) ...
            + 1/2*(3*lambda(p,3)-1).*(3*lambda(p,3)-2).*Dlambda3;
        w4(:,2*p-1:2*p) = 9/2*Dlambda3.*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*Dlambda2.*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambda2);
        w5(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,3).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*Dlambda3.*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambda3);
        w6(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,2).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambda2.*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambda1);
        w7(:,2*p-1:2*p) = 9/2*Dlambda3.*lambda(p,2).*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*Dlambda2.*(3*lambda(p,3)-1) ...
            + 9/2*lambda(p,3).*lambda(p,2).*(3*Dlambda3);
        w8(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,3).*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*Dlambda3.*(3*lambda(p,1)-1) ...
            + 9/2*lambda(p,1).*lambda(p,3).*(3*Dlambda1);
        w9(:,2*p-1:2*p) = 9/2*Dlambda1.*lambda(p,2).*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*Dlambda2.*(3*lambda(p,2)-1) ...
            + 9/2*lambda(p,1).*lambda(p,2).*(3*Dlambda2);
        w10(:,2*p-1:2*p) = 27*Dlambda1.*lambda(p,2).*lambda(p,3) ...
            + 27*lambda(p,1).*Dlambda2.*lambda(p,3) ...
            + 27*lambda(p,1).*lambda(p,2).*Dlambda3;
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dxx
if strcmpi(wStr,'.dxx')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = 9*Dlambdax(:,1).*Dlambdax(:,1).*(3*lambda(p,1)-1);
        w2(:,p) = 9*Dlambdax(:,2).*Dlambdax(:,2).*(3*lambda(p,2)-1);
        w3(:,p) = 9*Dlambdax(:,3).*Dlambdax(:,3).*(3*lambda(p,3)-1);
        w4(:,p) = 9/2*(6*Dlambdax(:,2).*Dlambdax(:,2)*lambda(p,3) + 6*Dlambdax(:,2).*lambda(p,2).*Dlambdax(:,3) ...
            + 6*lambda(p,2)*Dlambdax(:,2).*Dlambdax(:,3) - Dlambdax(:,2).*Dlambdax(:,3) - Dlambdax(:,2).*Dlambdax(:,3));
        w5(:,p) = 9/2*(6*Dlambdax(:,3).*Dlambdax(:,3)*lambda(p,1) + 6*Dlambdax(:,3).*lambda(p,3).*Dlambdax(:,1) ...
            + 6*lambda(p,3)*Dlambdax(:,3).*Dlambdax(:,1) - Dlambdax(:,3).*Dlambdax(:,1) - Dlambdax(:,3).*Dlambdax(:,1));
        w6(:,p) = 9/2*(6*Dlambdax(:,1).*Dlambdax(:,1)*lambda(p,2) + 6*Dlambdax(:,1).*lambda(p,1).*Dlambdax(:,2) ...
            + 6*lambda(p,1)*Dlambdax(:,1).*Dlambdax(:,2) - Dlambdax(:,1).*Dlambdax(:,2) - Dlambdax(:,1).*Dlambdax(:,2));
        w7(:,p) = 9/2*((6*lambda(p,3)-1)*(Dlambdax(:,2).*Dlambdax(:,3)+Dlambdax(:,3).*Dlambdax(:,2)) + 6*lambda(p,2)*Dlambdax(:,3).*Dlambdax(:,3));
        w8(:,p) = 9/2*((6*lambda(p,1)-1)*(Dlambdax(:,3).*Dlambdax(:,1)+Dlambdax(:,1).*Dlambdax(:,3)) + 6*lambda(p,3)*Dlambdax(:,1).*Dlambdax(:,1));
        w9(:,p) = 9/2*((6*lambda(p,2)-1)*(Dlambdax(:,1).*Dlambdax(:,2)+Dlambdax(:,2).*Dlambdax(:,1)) + 6*lambda(p,1)*Dlambdax(:,2).*Dlambdax(:,2));
        w10(:,p) = 27*(Dlambdax(:,1).*Dlambdax(:,2).*lambda(p,3) + Dlambdax(:,1)*lambda(p,2).*Dlambdax(:,3) + Dlambdax(:,1).*Dlambdax(:,2)*lambda(p,3) ...
            + lambda(p,1)*Dlambdax(:,2).*Dlambdax(:,3) + Dlambdax(:,1)*lambda(p,2).*Dlambdax(:,3) + lambda(p,1)*Dlambdax(:,2).*Dlambdax(:,3));
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dyy
if strcmpi(wStr,'.dyy')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = 9*Dlambday(:,1).*Dlambday(:,1).*(3*lambda(p,1)-1);
        w2(:,p) = 9*Dlambday(:,2).*Dlambday(:,2).*(3*lambda(p,2)-1);
        w3(:,p) = 9*Dlambday(:,3).*Dlambday(:,3).*(3*lambda(p,3)-1);
        w4(:,p) = 9/2*(6*Dlambday(:,2).*Dlambday(:,2)*lambda(p,3) + 6*Dlambday(:,2).*lambda(p,2).*Dlambday(:,3) ...
            + 6*lambda(p,2)*Dlambday(:,2).*Dlambday(:,3) - Dlambday(:,2).*Dlambday(:,3) - Dlambday(:,2).*Dlambday(:,3));
        w5(:,p) = 9/2*(6*Dlambday(:,3).*Dlambday(:,3)*lambda(p,1) + 6*Dlambday(:,3).*lambda(p,3).*Dlambday(:,1) ...
            + 6*lambda(p,3)*Dlambday(:,3).*Dlambday(:,1) - Dlambday(:,3).*Dlambday(:,1) - Dlambday(:,3).*Dlambday(:,1));
        w6(:,p) = 9/2*(6*Dlambday(:,1).*Dlambday(:,1)*lambda(p,2) + 6*Dlambday(:,1).*lambda(p,1).*Dlambday(:,2) ...
            + 6*lambda(p,1)*Dlambday(:,1).*Dlambday(:,2) - Dlambday(:,1).*Dlambday(:,2) - Dlambday(:,1).*Dlambday(:,2));
        w7(:,p) = 9/2*((6*lambda(p,3)-1)*(Dlambday(:,2).*Dlambday(:,3)+Dlambday(:,3).*Dlambday(:,2)) + 6*lambda(p,2)*Dlambday(:,3).*Dlambday(:,3));
        w8(:,p) = 9/2*((6*lambda(p,1)-1)*(Dlambday(:,3).*Dlambday(:,1)+Dlambday(:,1).*Dlambday(:,3)) + 6*lambda(p,3)*Dlambday(:,1).*Dlambday(:,1));
        w9(:,p) = 9/2*((6*lambda(p,2)-1)*(Dlambday(:,1).*Dlambday(:,2)+Dlambday(:,2).*Dlambday(:,1)) + 6*lambda(p,1)*Dlambday(:,2).*Dlambday(:,2));
        w10(:,p) = 27*(Dlambday(:,1).*Dlambday(:,2).*lambda(p,3) + Dlambday(:,1)*lambda(p,2).*Dlambday(:,3) + Dlambday(:,1).*Dlambday(:,2)*lambda(p,3) ...
            + lambda(p,1)*Dlambday(:,2).*Dlambday(:,3) + Dlambday(:,1)*lambda(p,2).*Dlambday(:,3) + lambda(p,1)*Dlambday(:,2).*Dlambday(:,3));
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end

%% u.dxy
if strcmpi(wStr,'.dxy') || strcmpi(wStr,'.dyx')
    [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10] = deal(zeros(NT,ng));
    for p = 1:ng
        w1(:,p) = 9*Dlambdax(:,1).*Dlambday(:,1).*(3*lambda(p,1)-1);
        w2(:,p) = 9*Dlambdax(:,2).*Dlambday(:,2).*(3*lambda(p,2)-1);
        w3(:,p) = 9*Dlambdax(:,3).*Dlambday(:,3).*(3*lambda(p,3)-1);
        w4(:,p) = 9/2*(6*Dlambdax(:,2).*Dlambday(:,2)*lambda(p,3) + 6*Dlambdax(:,2).*lambda(p,2).*Dlambday(:,3) ...
            + 6*lambda(p,2)*Dlambday(:,2).*Dlambdax(:,3) - Dlambdax(:,2).*Dlambday(:,3) - Dlambday(:,2).*Dlambdax(:,3));
        w5(:,p) = 9/2*(6*Dlambdax(:,3).*Dlambday(:,3)*lambda(p,1) + 6*Dlambdax(:,3).*lambda(p,3).*Dlambday(:,1) ...
            + 6*lambda(p,3)*Dlambday(:,3).*Dlambdax(:,1) - Dlambdax(:,3).*Dlambday(:,1) - Dlambday(:,3).*Dlambdax(:,1));
        w6(:,p) = 9/2*(6*Dlambdax(:,1).*Dlambday(:,1)*lambda(p,2) + 6*Dlambdax(:,1).*lambda(p,1).*Dlambday(:,2) ...
            + 6*lambda(p,1)*Dlambday(:,1).*Dlambdax(:,2) - Dlambdax(:,1).*Dlambday(:,2) - Dlambday(:,1).*Dlambdax(:,2));
        w7(:,p) = 9/2*((6*lambda(p,3)-1)*(Dlambdax(:,2).*Dlambday(:,3)+Dlambdax(:,3).*Dlambday(:,2)) + 6*lambda(p,2)*Dlambday(:,3).*Dlambdax(:,3));
        w8(:,p) = 9/2*((6*lambda(p,1)-1)*(Dlambdax(:,3).*Dlambday(:,1)+Dlambdax(:,1).*Dlambday(:,3)) + 6*lambda(p,3)*Dlambday(:,1).*Dlambdax(:,1));
        w9(:,p) = 9/2*((6*lambda(p,2)-1)*(Dlambdax(:,1).*Dlambday(:,2)+Dlambdax(:,2).*Dlambday(:,1)) + 6*lambda(p,1)*Dlambday(:,2).*Dlambdax(:,2));
        w10(:,p) = 27*(Dlambdax(:,1).*Dlambday(:,2).*lambda(p,3) + Dlambdax(:,1)*lambda(p,2).*Dlambday(:,3) + Dlambday(:,1).*Dlambdax(:,2)*lambda(p,3) ...
            + lambda(p,1)*Dlambdax(:,2).*Dlambday(:,3) + Dlambday(:,1)*lambda(p,2).*Dlambdax(:,3) + lambda(p,1)*Dlambday(:,2).*Dlambdax(:,3));
    end
    w = {w1,w2,w3,w4,w5,w6,w7,w8,w9,w10};
    return;
end