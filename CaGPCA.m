function [U1,V1] = CaGPCA(data, opts)
%------------------------------------------
% Date created: 03-01-2024
% @ Northwestern Polytechnical University 
% Please contact Lei Du and Jin Zhang(jinzhang@mail.nwpu.edu.cn) for any comments or questions.
% -----------------------------------------
X = data.X;
Y{1} = data.Y{1};
Y{2} = data.Y{2};
Y{3} = data.Y{3};
Y{4} = data.Y{4};
Y{5} = data.Y{5};
z = data.z;
n = size(X,1);
p = size(X,2);
q1 = size(Y{1},2);
q2 = size(Y{2},2);
q3 = size(Y{3},2);
q4 = size(Y{4},2);
q5 = size(Y{5},2);
U1 = ones(p, 5);
V1{1} = ones(q1, 1);
V1{2} = ones(q2, 1);
V1{3} = ones(q3, 1);
V1{4} = ones(q4, 1);
V1{5} = ones(q5, 1);
%% paramters initialization
W = ones(n,1)./n;
lambda0 = 0.1;
lambda1 = 0.1;
lambda2 = 0.1;
lambda5 = 0.1;
lambda_W = 0.001;
% set parameters
lambda_u1 = opts.lambda.u1;
lambda_u2 = opts.lambda.u2;
lambda_u3 = opts.lambda.u3;
lambda_v1 = opts.lambda.v1;
lambda_v2 = opts.lambda.v2;
lambda_v3 = opts.lambda.v3;
lambda_v4 = opts.lambda.v4;
lambda_v5 = opts.lambda.v5;
% set stopping criteria
max_Iter = 100;
t = 0;
tol = 1e-5;
tu1 = inf;
tv1 = inf;
% initial scale
for i = 1 : 5
    U = U1(:, i) ;
    scale = sqrt(U' * X'*X * U);
    U1(:, i) = U1(:, i) ./ scale;
    V = V1{i};
    scale = sqrt(V'* Y{i}'* Y{i}* V);
    V1{i} = V1{i} ./ scale;
end
W = ones(n,1);
block = 50;
nblock = ceil(p/block);
Xall = X;
Yv1 = Y{1}*V1{1};
Yv2 = Y{2}*V1{2};
Yv3 = Y{3}*V1{3};
Yv4 = Y{4}*V1{4};
Yv5 = Y{5}*V1{5};
while (t<max_Iter && (tu1>tol|| tv1>tol))
    t = t + 1;
    U1_old = U1;
    V1_old = V1;
    Yv = [Yv1 Yv2 Yv3 Yv4 Yv5];
        ut = [];
        su1 = 0;
        su2 = 0;
        su3 = 0;
        su4 = 0;
        su5 = 0;
        U_old = U1;       
        for iu = 1:nblock
        if iu*block <= p
            X = Xall(:,1+(iu-1)*block:iu*block);
            sub_u = U1(1+(iu-1)*block:iu*block,:);
        else
            X = Xall(:,1+(iu-1)*block:end);
            sub_u = U1(1+(iu-1)*block:end,:);
        end
        XX = X'*diag(W) * X;
        for i = 1 : 5
            XY{i} = X' *diag(W)* Y{i};
            YY{i} = Y{i}'*diag(W)* Y{i};
            YX{i} = XY{i}';
        end
        d1 = updateD(sub_u);
        D1 = diag(d1);
        fd1 = updateD_FGL21(sub_u(:, 1), sub_u(:, 2), sub_u(:, 3), sub_u(:, 4),sub_u(:, 5));
        fD1 = diag(fd1); 
        for i = 1 : 5
            du1 = updateD(sub_u(:,i));
            Du1{i} = diag(du1);
        end
        for i = 1 : 5
        F1 = XX + lambda_u1 * D1 + lambda_u2 *fD1 + lambda_u3 * Du1{i};
        b1 = XY{i} * (V1{i});
        sub_u(:,i) = F1 \ b1;
        end
        ut = [ut; sub_u];
        [su1,su2,su3,su4,su5] = cal_scale(sub_u,XX,su1,su2,su3,su4,su5);
    end
     U1(:, 1) = ut(:,1) / sqrt(su1);
     U1(:, 2) = ut(:,2) / sqrt(su2);
     U1(:, 3) = ut(:,3) / sqrt(su3);
     U1(:, 4) = ut(:,4) / sqrt(su4);
     U1(:, 5) = ut(:,5) / sqrt(su5);
    for i = 1 : 5
        XY{i} = Xall' *diag(W)* Y{i};
        YY{i} = Y{i}'*diag(W)* Y{i};
        YX{i} = XY{i}';
    end
        d21 = updateD(V1{1});
        D21 = diag(d21);
        F1 = YY{1} + lambda_v1 * D21;
        b1 = YX{1} * U1(:, 1) + Y{1}'*z;
        V1{1} = F1 \ b1;
        V = V1{1};
        scale = sqrt(V' * Y{1}' * Y{1} * V);
        V1{1} = V1{1} ./ scale;
        d22 = updateD(V1{2});
        D22 = diag(d22);
        F2 = YY{2} + lambda_v2 * D22;
        b2 = YX{2} * U1(:, 2) + Y{2}'*z;
        V1{2} = F2 \ b2;
        V = V1{2};
        scale = sqrt(V' * Y{2}' * Y{2} * V);
        V1{2} = V1{2} ./ scale;
        d23 = updateD(V1{3});
        D23 = diag(d23);
        F2 = YY{3} + lambda_v3 * D23;
        b2 = YX{3} * U1(:, 3) + Y{3}'*z;
        V1{3} = F2 \ b2;
        V = V1{3};
        scale = sqrt(V' * Y{3}' * Y{3} * V);
        V1{3} = V1{3} ./ scale;
        
        d24 = updateD(V1{4});
        D24 = diag(d24);
        F2 = YY{4} + lambda_v4 * D24;
        b2 = YX{4} * U1(:, 4 ) + Y{4}'*z;
        V1{4} = F2 \ b2;
        V = V1{4};
        scale = sqrt(V' * Y{4}' * Y{4} * V);
        V1{4} = V1{4} ./ scale;
        d25 = updateD(V1{5});
        D25 = diag(d25);
        F2 = YY{5} + lambda_v5 * D25;
        b2 = YX{5} * U1(:, 5) + Y{5}'*z;
        V1{5} = F2 \ b2;
        V = V1{5};
        scale = sqrt(V' * Y{5}' * Y{5} * V);
        V1{5} = V1{5} ./ scale;
        
    obj = norm(Xall*U1(:,1)-Y{1}*V1{1})^ 2 + norm(Xall*U1(:,2)-Y{2}*V1{2})^ 2 + norm(Xall*U1(:,3)-Y{3}*V1{3})^ 2 ... 
    + norm(Xall*U1(:,4)-Y{4}*V1{4})^ 2 + norm(Xall*U1(:,5)-Y{5}*V1{5})^ 2;
    grad_W = 2 * lambda0 * obj.*W...
            + lambda1 * balance_gradfast(W, Xall) * ones(p,1)...
            + 4 * lambda2*W.*W.*W...           
            + 4*lambda5*(sum(W.*W)-1)*W;
        W = W - lambda_W*grad_W;
        W = W./sqrt(sum(W.*W));
        W = W.*W;  
        
    tu1 = max(max(abs(U1 - U1_old)));
    t1 = max(max(abs(V1{1} - V1_old{1})));
    t2 = max(max(abs(V1{2} - V1_old{2})));
    t3 = max(max(abs(V1{3} - V1_old{3})));
    t4 = max(max(abs(V1{3} - V1_old{3})));
    t5 = max(max(abs(V1{3} - V1_old{3})));
    t6 = max(max(abs(V1{3} - V1_old{3})));
    tv1 =max([t1,t2,t3,t4,t5,t6]);
     end
end

function [sum1 sum2 sum3 sum4 sum5] = cal_scale(beta,XX,sum1,sum2,sum3,sum4,sum5)
[p,ntask] = size(beta);
if nargin <3
    for i = 1:ntask
        temp1 = beta(:,i);
        temp2 = temp1';
        eval(['sum',str2double(i),' = temp2*XX*temp1;'])
    end
else
    for i = 1:ntask
        temp1 = beta(:,i);
        temp2 = temp1';
        eval(['sum',num2str(i),' = sum',num2str(i),'+temp2*XX*temp1;']);
    end
end
end

function D = updateD(W, group)
    [n_features, n_tasks] = size(W);
        for i = 1 : n_features
            d(i) = sqrt(sum(W(i, :) .^ 2) + eps);
        end
    D = 0.5 ./ d;
end

function [D] = updateD_FGL21(u1, u2, u3, u4, u5)
ulen = length(u1);
for i = 1:ulen
    if i == 1
        d(i) = sqrt(u1(i).^2+u2(i).^2+u3(i).^2+u4(i).^2+u5(i).^2+u1(i+1).^2+u2(i+1).^2+u3(i+1).^2+u4(i+1).^2+u5(i+1).^2+eps);
        d(i) = 0.5 ./ d(i);
    elseif i==ulen
        d(i) = sqrt(u1(i-1).^2+u2(i-1).^2+u3(i-1).^2+u4(i-1).^2+u5(i-1).^2+u1(i).^2+u2(i).^2+u3(i).^2+u4(i).^2+u5(i).^2+eps);
        d(i) = 0.5 ./ d(i);
    else
        d(i) = 0.5./(sqrt(u1(i-1).^2+u2(i-1).^2+u3(i-1).^2+u4(i-1).^2+u5(i-1).^2+u1(i).^2+u2(i).^2+u3(i).^2+u4(i).^2+u5(i).^2+eps))+0.5./(sqrt(u1(i).^2+u2(i).^2+u3(i).^2+u4(i).^2+u5(i).^2+u1(i+1).^2+u2(i+1).^2+u3(i+1).^2+u4(i+1).^2+u5(i+1).^2+eps));
    end
    D = d;
end
end




