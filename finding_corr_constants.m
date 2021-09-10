function[A,a_ij,b_ij]=finding_corr_constants(X,U,a,b,c,d,noise_var,lambda)

C = X*X'/length(X)-noise_var*eye(size(X,1));

U = kron(U,U);

a_ij = -((a*a')./(c*ones(1,length(c))+ones(length(c),1)*c'));
a_ij = a_ij - ((a*b')./(c*ones(1,length(c))+ones(length(c),1)*d'));

b_ij = -((b*a')./(d*ones(1,length(c))+ones(length(c),1)*c'));
b_ij = b_ij - ((b*b')./(d*ones(1,length(d))+ones(length(d),1)*d'));

K0 = tril(a_ij + b_ij) + triu(a_ij' + b_ij')-diag(diag(a_ij + b_ij));
K0 = K0(:);

select = sqrt(lambda)*ones(1,length(lambda))-ones(length(lambda),1)*sqrt(lambda');
select = exp(-2*(select.^2));
select = select(:);
select = select .* (select > 10^(-5));

[ind,~] = find(select);

ind_diag = 1:length(c)+1:(length(c))^2;
ind = unique([ind;ind_diag']);

[diag_entries,~] = find(ind==ind_diag);
A = sdpvar(length(c),length(c));
a_vec = A(:);

SA = eye(length(c)^2);
SA = SA(:,ind);

Obj   = (norm(U*SA*diag(K0(ind))*a_vec(ind)-C(:)))^2;
Const = (A >= 0);
ops = sdpsettings('solver','sedumi','debug',1);
sol = optimize(Const,Obj,ops)

A = value(A);

% B = zeros(length(c)^2,1);
% % A = A.*select(ind);
% B(ind) = A;
% 
% A =B;
% 
% A = reshape(A',length(c),length(c));
% 
% 
A = A.*(abs(A)>(10^(-8)));
A(isnan(A))=0;
A(isinf(A))=0;

if sol.problem==1
   A=zeros(size(A)); 
end
end