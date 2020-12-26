clearvars, clearvars –global 
clear all, clear global 

% ======= variables ========
global dim R m P_0 basis;
dim = 25;                                        
R = 10;             
m = 2;                 
P_0 = 339.3 ;  

basis = cell(dim, 3);

for n = 1:1:dim
    basis{n, 1} = chebfun(@(r) sin(r * n * pi ./ R ), [0, R]);
    gram_schmidt(n);  
    basis{n, 2} = diff(basis{n, 1});
    basis{n, 3} = diff(basis{n, 2});
end

t = cputime;
disp(['dim: ',num2str(dim), ...
      ' R: ',num2str(R),    ...
      ' m: ',num2str(m),    ...
      ' P_0: ',num2str(P_0)    ...
      ]) 

% ============ problem setup ================
c = optimvar('c', dim);
obj = fcn2optimexpr(@I, c); 
prob = optimproblem('Objective', obj);
cons1 = sum(c .* c) == P_0;          % beam power constraints

prob.Constraints.cons1 = cons1; 

options = optimoptions(@fmincon,...
                        'Algorithm','active-set',...
                        'MaxFunEvals', 1e+05, ...
                        'MaxIter', 1000,...
                        'TolCon',1e-8,...
                        'TolFun',1e-8);

x0.c =  zeros(1, dim);
x0.c(1) = 1;

try
    [sol,fval,exitflag,output] = solve(prob,x0,'options',options);
    % display success/fail
    disp(exitflag > 0)
    if exitflag > 0 % if successful optimization

        % plotting 
        u = lin_comb(sol.c);
        cache = sol.c;

        % compute eigenvalue
        eigen = compute_eigen(sol.c, P_0);
        disp(['Yields a good eigenvalue:', num2str(eigen)])

        % compute error 
        error = compute_error(sol.c, eigen);
        disp(['error:', num2str(error)])
    end
catch
    disp('some error')
end
disp(['Time it took:', num2str(cputime-t)])


% ================== functions ===================

% Linear Combination of variational vector and basis (global)
function u = lin_comb(vec)
    global R dim basis;
    u = chebfun(0, [0 R]);
    for n = 1:dim
        u = u + vec(n) * basis{n, 1};
    end
end


function u_r = lin_comb_r(vec)
    global R dim basis;
    u_r = chebfun(0, [0 R]);
    for n = 1:dim
        u_r = u_r + vec(n) * basis{n, 2};
    end
end

function u_rr = lin_comb_rr(vec)
    global R dim basis;
    u_rr = chebfun(0, [0 R]);
    for n = 1:dim
        u_rr = u_rr + vec(n) * basis{n, 3};
    end
end


% inner product defined as: 2 pi * int_0^R r*f(r)*g(r) dr
function prod = inner_product(func1, func2)
    global R;
    prod = func1 * func2 * chebfun(@(r) r, [0, R]);
    prod = 2 * pi * sum(prod);
end

% gram_schmidt process on n th basis
function gram_schmidt(n)
    global basis;
    basis{n, 1} =  basis{n, 1} / sqrt(inner_product(basis{n, 1}, basis{n, 1}));
    for i = 1 : n-1
        basis{n, 1} = basis{n, 1} - projection(basis{i, 1}, basis{n, 1});
    end
    basis{n, 1} = basis{n, 1} / sqrt(inner_product(basis{n, 1}, basis{n, 1}));
    
    function proj = projection(fun_i, fun_n)
        proj = inner_product(fun_i,fun_n) * basis{i, 1} / inner_product(fun_i,fun_i);
    end

end

% I: action functional
function integral = I(vec)
    global R m;
    u = lin_comb(vec);
    u_r = lin_comb_r(vec);
    r = chebfun(@(r) r, [0, R]);
    integrand = r * u_r^2;
    integrand = integrand + m^2 * u^2 / r;
    integrand = integrand + 4 * r * sqrt(1 + u^2);
    integral = sum(integrand) / 2;
end 

function eigen = compute_eigen(vec, p_0)
    global R m;
    u = lin_comb(vec);
    u_r = lin_comb_r(vec);
    r = chebfun(@(r) r, [0, R]);
    ut = r * u_r^2;
    ut = ut + m^2 * u^2 / r;
    ut = ut + 2 * r * u^2 / sqrt(1 + u^2);
    eigen = -2 * pi * sum(ut)/ p_0;
end

% compute error intergrating (1)
function error = compute_error(vec, b)
    global R m;
    u = lin_comb(vec);
    u_r = lin_comb_r(vec);
    r = chebfun(@(r) r, [0, R]);
    ut = 2 * b * r * u;
    ut = ut - u_r - r * lin_comb_rr(vec);
    ut = ut + m^2 * u / r ;
    ut = ut - 2 * r * u ;
    ut = ut + 2 * r * u / sqrt(1 + u^2);
    error = sum(ut^2);
end
