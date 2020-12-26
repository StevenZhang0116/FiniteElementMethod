filename = 'eigenrmp6.mat'; 
eigen_mat = matfile(filename).eigen_mat;
% -0.50966
% 0.00056732
for i =1:40
    disp(i)
    disp(eigen_mat(1, 1, i, 1))
end

save(filename, 'eigen_mat')