% This code calculates a distribution of relaxation times (DRT) for a given
% EIS data.

clear; clc; close all
%% Configuration

% Data

    % Data path and name
    %filename_data = "1. Li'SE'Li after 0.5 mAh plating GEIS.csv";
    %filename_data = "1. Li'SE'Li after 0.5 mAh plating PEIS.csv";
    filename_data = "2. SUS'SE'Li after 0.5 mAh plating PEIS.csv";
    
   % Freq range to include
    f_range = [-inf inf]; % [-inf inf] to include all data

   
%  Config DRT

   % Hyperparameters
    l_star = -4; % regularization parameter (lambda), log scale
                % --> 적절한 l (lambda) 값을 구하기 위해서는 run_DRT_CV.m 파일 이용.

    t_lb = -6; % time scale to include in DRT analysis, lower limit, log scale   
    t_ub = 4; % time scale to include in DRT analysis, upper limit, log scale
    ppd_t = 10; % point per decade in t-domain



    % Hard setting (바꾸기 전에 문의 요 - 송주현)
        % Treating the intercept resistance 
        type_intercept = 1; % 0, use given R_intercept / 1, include fitting R_intercept in the DRT inversion
        % Type of kernel function
        type_kernel = 0; % 0,RC (DRT) / 1,FL (DDT) / 2,BD  (DDT)
                % 0 번을 지정하면 DRT 로 작동합니다 - 현재 코드는 DDT 비활성화됨.
        % Type of weighting 
        type_weight = 0; % 0, Uniform / 1, Relative 


%% Data Pre-processing

    % Load the selected data
    data = readmatrix(filename_data); % 숫자만 읽어 옴. expecting [freq (Hz), Z_real (Ohm), Z_ohm]

    % 지정된 주파수 범위 외 삭제
    data = data(data(:,1) > f_range(1),:);
    data = data(data(:,1) < f_range(2),:);

    % 허수 임피던스 음수 부분 (고주파수 인덕턴스) 빼기
    data = data(data(:,3)>0,:);

    % DRT 변환에 포함될 데이터 지정하기
    f_data = data(:,1); % [Hz] frequency data (vector)
    z_data = data(:,2) -data(:,3)*1i; % [Ohm] impedance data (vector, complex)
    


    % plot the imported impedance data
    figure(1);
    subplot(1,3,1)
    hold on;
    plot(real(z_data),-imag(z_data),'ko')
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0.8*min(real(z_data)), 1.2*max(real(z_data))],'Ylim',[0,1.2*max(real(z_data))-0.8*min(real(z_data))])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    legend_fig1{1} = 'exp data';
    legend(legend_fig1,'location','northwest')



%% Preparing the Inversion

    % Define the t and tau 
       % generate t vector and tau vector
        t = t_lb:1/ppd_t:t_ub;
        tau = exp(t); % tau is define in Reference page 2, t = log(tau)
        N = length(t);

    % Calculate the Kernels
        % frequency and x
        w = 2*pi*f_data;
        x = -log(w); % x is defined in Reference page 2, x = -log(w)
        M = length(x);

        % Kernel
        [K,~,D] = DRT_Kernel(type_intercept,type_kernel,x,ppd_t,t,zeros(M,1));




%% DRT Inversion: (3) Inversion with the optimal hyperparameters

    % Calculate weighting matrix
        [W_re,W_im] = DRT_Weight(z_data,M,type_weight);
    % Inversion by quadratic programming
        q_star = DRT_Inversion(z_data,t,K,W_re,W_im,D,l_star,1,1,type_kernel);
    % Fitted model evaluation
        z_D_star = K*q_star;
    % Normalization
       % q_star = q_star(2:end)./sum(q_star(2:end))*ppd_t;


    % Plot the solution
        figure(1); subplot(1,3,1); hold on;
        plot(real(z_D_star),-imag(z_D_star),'-b','linewidth',2);
        legend_fig1{end+1} = 'DRT model';
        legend(legend_fig1)

        figure(1); subplot(1,3,2); hold on;
        plot(t,q_star(2:end),'-b')
        axis([t_lb,t_ub,0,1.1*max(q_star(2:end-1))]);
        xlabel ('t = log(\tau)')
        ylabel ('q = \tauP(\tau)')
        legend_fig4{1} = 'DRT solution';
        legend(legend_fig4)

    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')

    set(gcf,'position',[100 600 900 400])

% end.