clearvars, clearvars –global 
clear all, clear global 



% ========= copy to plot ========
filename = 'eigenrmp6.mat';         
dim = 30;                                         
R_range = [20];                   
m_range = [1,2,3,4,5];                 
P_0_range = [1,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000]               
% %===============================

% nr = length(R_range);
% nm = length(m_range);
% np = length(P_0_range);
% 
% eigen_mat = matfile(filename).eigen_mat;
% 
% r_range = 1:nr;
% p_range = 1:np;
% disp(eigen_mat)
% 
% for r = r_range
%     loglog(P_0_range, reshape(eigen_mat(r, 1, :, 1), 1, [])), hold on
% end
% ylabel('Beta, ')
% xlabel(['P_{0}'])

% hold off
% for p_0 = P_0_range
%     plot(R_range, reshape(eigen_mat(:, 1, p_0), 1, [])), hold on
%     plot(R_range, reshape(eigen_mat(:, 1, p_0), 1, [])), hold on
% end
% 
% ylabel('Beta')
% xlabel('R')
% legend(string(P_0))
% plot(R, reshape(eigen_mat(:, 1, 1), 1, []), 'b')
% 
% plot(P_0_range, reshape(eigen_mat(3, 1, :, 1), 1, [])), hold on
% 
% 
% for P = p_range
%     errorbar(R_range,...
%              reshape(eigen_mat(:, 1, P, 1), 1, []), ...
%              reshape(eigen_mat(:, 1, P, 2), 1, []), '-s')
%     hold on
% end
% legend(string(P_0_range))
% ylabel('Beta')
% xlabel('R')
% ylim([-20, 0])
% xlim([0, 20])
% 
% for m = m_range
% %     plot(P_0_range, reshape(eigen_mat(1, m, :, 1), 1, []), '-s')
%     errorbar(P_0_range,...
%              reshape(eigen_mat(1, m, :, 1), 1, []), ...
%              reshape(eigen_mat(1, m, :, 2), 1, []), '-s')
%     hold on
%     
% %     plot(P_0_range, ones(size(P_0_range)) * -m^2/400 - 1)
% %     hold on
% end
% ylim([-2 0])
% xlim([0 11000])
% legend('m=1', 'm=2', 'm=3', 'm=4', 'm=5')
% ylabel('\beta, Propagation Constant')
% xlabel('P_{0}, Beam Power functional')
% 
% for m = m_range
%     errorbar(R_range,...
%              reshape(eigen_mat(:, m, 1, 1), 1, []), ...
%              reshape(eigen_mat(:, m, 1, 2), 1, []), '-s')
%     hold on
% end
% 
% legend(strcat('m=',string(m_range)))
% ylabel('\beta, Propagation Constant')
% xlabel('Radius')
% ylim([-2, 0])
% xlim([-5, 3.5])
