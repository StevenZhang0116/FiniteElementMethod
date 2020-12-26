clearvars, clearvars –global 
clear all, clear global 

%  parallel computign setup 
delete(gcp('nocreate'))
parpool('local')
UseParallel = true;

% ========== variables ===========
global dim R m P_0 basis;

% copy and paste this box to "plottingpad.m"
% ========= copy to plot ========
filename = 'eigenrmp5.mat';          
dim = 30;                   
R_range = [20];                 
m_range = [1,2,3,4,5,6,7,8];                 
P_0_range = [150];                
% %===============================


basis = cell(dim, 3);
nr = length(R_range);
nm = length(m_range);
np = length(P_0_range);

%{
eigen values stored in a rank 4 tensor with 
axis 1: index for R
axis 2: index for m
axis 3: index for P_0
axis 4: value and corresponding error and solution vector. 
        with 1 being eigenvalue, and 
        2 being its corresponding error, and
        3 to 2 + dim being the solution to variational vector
%}


eigen_mat = NaN(nr, nm, np, 2 + dim);

%{
cache is used to set initial point in fmincon
solution of variational vector of previous iteration
is set to to be initial point. set to (x01,0,0,...) if
r is changed
%}
cache = zeros(1, dim);
x01 = 1;

% index used for for-loop. Can be used to pick up half finished code

m_ind = 1;
p_ind = 1;
r_ind = 1;


for R = R_range(r_ind:end)
    
    % ======= creating basis ===========
    for n = 1:1:dim
        basis{n, 1} = chebfun(@(r) sin(r * n * pi ./ R ), [0, R]);
        gram_schmidt(n);  
        basis{n, 2} = diff(basis{n, 1});
        basis{n, 3} = diff(basis{n, 2});
    end
    cache(1) = x01;  % initial value
    
    for m = m_range(m_ind:end)
        for P_0 = P_0_range(p_ind:end)
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
            x0.c = cache;

            % ============== optimization ==============='
            try
                [sol,fval,exitflag,output] = solve(prob, x0 ,'options', options);
                % display success/fail
                disp(exitflag > 0)
                if exitflag > 0 % if successful optimization

                    % plotting 
                    u = lin_comb(sol.c);
                    cache = sol.c;

                    % compute eigenvalue
                    eigen = compute_eigen(sol.c, P_0);
                    eigen_mat(r_ind, m_ind, p_ind, 1) = eigen;
                    disp(['Yields a good eigenvalue:', num2str(eigen)])

                    % compute error 
                    error = compute_error(sol.c, eigen);
                    eigen_mat(r_ind, m_ind, p_ind, 2) = error;
                    disp(['error:', num2str(error)])
                    
                    eigen_mat(r_ind, m_ind, p_ind, 3:2+dim) = sol.c;
                    save(filename, 'eigen_mat')
                end
            catch
                disp('some error')
            end
            disp(['Time it took:', num2str(cputime-t)])
            p_ind = p_ind + 1;
        end
        p_ind = 1;
        m_ind = m_ind + 1;
        cache = zeros(1, dim);
        cache(1) = x01;
    end
    r_ind = r_ind + 1;
    p_ind = 1;
    m_ind = 1;
end
save(filename, 'eigen_mat')


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
    integrand = integrand - 2 * r * u^2;
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
    ut = ut - 2 * r * u^2;
    eigen = -pi * sum(ut)/ p_0;
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

