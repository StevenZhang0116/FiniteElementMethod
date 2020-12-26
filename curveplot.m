clearvars, clearvars –global 
clear all, clear global 
% ======= variables ========
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

eigen_mat = matfile(filename).eigen_mat;

r_ind = 1;
m_ind = 1;
p_ind = 1;

for R = R_range(r_ind:end)
    
    % ======= creating basis ===========
    for n = 1:1:dim
        basis{n, 1} = chebfun(@(r) sin(r * n * pi ./ R ), [0, R]);
        gram_schmidt(n);  
        basis{n, 2} = diff(basis{n, 1});
        basis{n, 3} = diff(basis{n, 2});
    end

    for m = m_range(m_ind:end)
        for P_0 = P_0_range(p_ind:end)
            u = lin_comb(eigen_mat(r_ind, m_ind, p_ind, 3:2+dim));
            disp(eigen_mat(r_ind, m_ind, p_ind, 1))
            disp(eigen_mat(r_ind, m_ind, p_ind, 2))
            plot(u,'LineWidth',1.5), hold on
            
            p_ind = p_ind + 1;
        end
        p_ind = 1;
        m_ind = m_ind + 1;
        cache = zeros(1, dim);
    end
    r_ind = r_ind + 1;
    p_ind = 1;
    m_ind = 1;
end
ylabel('u(r), Soliton aplitude')
xlabel('r, Distance from vortex core')
xlim([0, 22])



% ================== functions ===================

% Linear Combination of variational vector and basis (global)
function u = lin_comb(vec)
    global R dim basis;
    u = chebfun(0, [0 R]);
    for n = 1:dim
        u = u + vec(n) * basis{n, 1};
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
