function F = benchmark_functions(node, Q,  Ri, A,x)
    xopt=x;

    id = node.id;
    switch node.funid
        case 1
            F = f1_Elliptic_shifted_rot(Q,  Ri,xopt,A,id);

        case 2
            F = f2_Schwefel_shifted_rot(Q,  Ri,xopt,A,id);

        case 3
            F = f3_Rosenbrock_shifted_rot(Q,  Ri,xopt,A,id);

        case 4
            F = f4_Griewank_shifted_rot(Q,  Ri,xopt,A,id);

        case 5
            F = f5_Elli_Schw_shifted_rot(Q,  Ri,  id,xopt);

        case 6
            F = f6_Elli_Rosen_shifted_rot(Q,  Ri,  id,xopt);

        case 7
            F = f7_Schw_Rosen_shifted_rot(Q,  Ri,  id,xopt);

        case 8
            F = f8_Elli_Griew_shifted_rot(Q,  Ri,  id,xopt);

        case 9
            F = f9_Schw_Griew_shifted_rot(Q,  Ri,  id,xopt);

        case 10
            F = f10_Rosen_Griew_shifted_rot(Q,  Ri,  id,xopt,A);
    end
end
%基础函数
function [f_elem, Z] = elem_elliptic(Q, xopt, Ri)
    [~, D] = size(Q);

    % Shift
    Q_shift = Q - xopt;

    % Rotation
    Z =  (Ri*Q_shift')';

    % Osz + Asy + Lambda
    Z = Tosz(Z);


    % Elliptic 权重: 10^6^(i/(D-1)) = 10.^((0:D-1)/(D-1)*6)
    w = 10.^((0:D-1)/(D-1)*6);     % 1×D
    f_elem = sum(Z.^2 .* w, 2);    % N×1
end

function [f_elem, Z] = elem_schwefel(Q, xopt, Ri)
  

    % Shift
    Q_shift = Q - xopt;

    % Rotation
    Z =  (Ri*Q_shift')';

    % Osz + Asy（无 Lambda）
    Z = Tosz(Z);
    
    Z = Tasy(Z, 0.2);

    % Schwefel 1.2: sum( (cumsum(Z,2)).^2 , 2 )
    f_elem = sum(cumsum(Z, 2).^2, 2);
end

function [f_elem, Z] = elem_rosenbrock(Q, xopt, Ri)
    [N, ~] = size(Q);

    % Shift
    Q_shift = Q - xopt;

    % Rotation
    Z =  (Ri*Q_shift')';

    % Rosenbrock: sum_{j=1}^{D-1} [100 (z_{j+1} - z_j^2)^2 + (z_j - 1)^2]
    f_elem = sum( 100*(Z(:,2:end) - Z(:,1:end-1).^2).^2 + ...
                  (Z(:,1:end-1) - 1).^2, 2 );
end

function [f_elem, Z] = elem_griewank(Q, xopt, Ri)
    [N, D] = size(Q);

    % Shift
    Q_shift = Q - xopt;

    % Rotation
    Z = Q_shift * Ri;
    
    % Griewank: sum(z^2)/4000 - prod(cos(z/sqrt(j))) + 1
    for node=1:N
    end
    sum_term = sum(Z.^2, 2) / 4000;
    cos_term = prod(cos(Z ./ sqrt(1:D)), 2);
    f_elem   = sum_term - cos_term + 1;
end

%同构函数
function f = f1_Elliptic_shifted_rot(Q,  Ri,xopt,A,id)
    [f_elem, Z] = elem_elliptic(Q, xopt, Ri);

    % f_linear=linear(Z,A,100,id);

    % f = f_elem + f_linear ;
    f = f_elem ;
end
function f = f2_Schwefel_shifted_rot(Q,  Ri,xopt,A,id)

    [f_elem, Z] = elem_schwefel(Q, xopt, Ri);
    
   

    f = f_elem ;
end
function f = f3_Rosenbrock_shifted_rot(Q,  Ri,xopt,A,id)

    [f_elem, Z] = elem_rosenbrock(Q, xopt, Ri);
      
    % f_linear=linear(Z,A,100,id);

    % f = f_elem + f_linear ;
    f=f_elem;
end
function f = f4_Griewank_shifted_rot(Q, Ri,xopt,A,id)

    
    [f_elem, Z] = elem_griewank(Q, xopt, Ri);
     
    

    f = f_elem  ;
end

%混合函数  mod 0 Elli  mod 1 Schwefel
function f = f5_Elli_Schw_shifted_rot(Q,  Ri,  id,xopt)

    idx0 = mod(id-1, 2);   % 0 -> 第1个函数, 1 -> 第2个函数

    if idx0 == 0
      
        [f_elem, ~] = elem_elliptic(Q, xopt, Ri);
    else

        [f_elem, ~] = elem_schwefel(Q, xopt, Ri);
    end

    
    f = f_elem;
end
function f = f6_Elli_Rosen_shifted_rot(Q,  Ri,  id,xopt)
 
    idx0 = mod(id-1, 2);

    if idx0 == 0
        [f_elem, ~] = elem_elliptic(Q, xopt, Ri);
    else
        [f_elem, ~] = elem_rosenbrock(Q, xopt, Ri);
    end

 

    f = f_elem ;
end
function f = f7_Schw_Rosen_shifted_rot(Q,  Ri,  id,xopt)
   
    idx0 = mod(id-1, 2);
  
    if idx0 == 0
        [f_elem, ~] = elem_schwefel(Q, xopt, Ri);
    else
        [f_elem, ~] = elem_rosenbrock(Q, xopt, Ri);
    end

 
    f = f_elem;
end
function f = f8_Elli_Griew_shifted_rot(Q,  Ri,  id,xopt)

    idx0 = mod(id-1, 2);

    if idx0 == 0

        [f_elem, ~] = elem_elliptic(Q, xopt, Ri);
    else

        [f_elem, ~] = elem_griewank(Q, xopt, Ri);
    end

   
    f = f_elem ;
end
function f = f9_Schw_Griew_shifted_rot(Q, Ri,  id,xopt)
    idx0 = mod(id-1, 2);
    if idx0 == 0
        [f_elem, ~] = elem_schwefel(Q, xopt, Ri);
    else
        [f_elem, ~] = elem_griewank(Q, xopt, Ri);
    end

   
    f = f_elem ;
end
function f = f10_Rosen_Griew_shifted_rot(Q,  Ri,  id,xopt,A)
 
    idx0 = mod(id-1, 2);
    if idx0 == 0
        [f_elem, Z] = elem_rosenbrock(Q, xopt, Ri);
    else
        [f_elem, Z] = elem_griewank(Q, xopt, Ri);
    end

    % f_linear=linear(Z,A,100,id);

    f = f_elem  ;
end

function f=linear(Z,A,Weight,id)
    
    % 假设 Z 是一个 n x m 的矩阵，A 是一个矩阵，id 是行索引，weight 是加权因子
    f = Weight * sum(Z .* repmat(A(id, :), size(Z, 1), 1), 2);


end

% 变换函数
function Y = Tosz(X)
    [N,M] = size(X);
    for i = 1:N
        for j = 1:M
            t = X(i,j);
            
            if t>0
                 c1=10;c2=7.9;
            else 
                 c1=5.5;c2=3.1;
            end
            y = log(abs(t)) + 0.049 * ( ...
                sin(c1 * log(abs(t))) + ...
                sin(c2 * log(abs(t))) );
            
            X(i,j) = sign(t) * exp(y);
        end
    end
    Y = X;
end

function Y = Tasy(X, beta)
    if nargin < 2
        beta = 0.2;
    end
    [N,M] = size(X);
    for i = 1:N
        for j = 1:M
            t = X(i,j);
            if t > 0
                X(i,j) = t^(1 + beta * (j-1)/(M-1) * sqrt(t)); % 注意：索引从1开始，所以用 j-1
            end
        end
    end
    Y = X;
end

function Y = LambdaX(X, alpha)
    % 对应 C++ Lambda(z, alpha, dim):
    % z[i] = z[i] * alpha^(0.5 * i/(dim-1))
    if nargin < 2
        alpha = 10;
    end
    [N,M] = size(X);
    exponents = 0.5 * (0:M-1) / (M-1);      % 1×M
    scales = alpha .^ exponents;            % 1×M
    Y = X .* scales;                        % N×M，每列按比例缩放
end
