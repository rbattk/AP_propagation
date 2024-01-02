%Biyolektrik ve Biyomanyetizma Dersi Proje Ödevi                          %
%Aksiyon Potansiyelinin Oluşumu ve Yayılımı 2.Aşama                       %
%Rabia TUTUK ---> 36203630007                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
clear all;
 
V_m = 0; V_r = -90; V_Na = 115; V_K = -12; V_L = 10.613; 
G_Na_max = 120; G_K_max = 36; G_L_max = 0.3; 
ri = 30000;ro = 1000; a = 0.01;
c_m = 1; C_m = 2*pi*(ri+ro)*a*c_m;
a_p_p_t = 40; aralik=0.01;t=0:aralik:a_p_p_t; k=numel(t);
u_ar = 1; L=10; l = 0:u_ar:L; p = numel(l);

%I(1:k) = [100]; 
%f = iki uyaran arasındaki süre
f = input('İki uyaran arasındaki zaman farkını giriniz:');
If = input('Birinci uyaranın genliğini giriniz:');
I = input('İkinci uyaranın genliğini giriniz:');

f = 1:1:f;
I(1:100) = If; I(101:numel(f)) = 0; I(numel(f)+1 : numel(f)+100)=If; I(numel(f)+101:k)=0;

alpha_n = (0.1 - (0.001 * V_m)) / (exp(1 - (0.1 * V_m)) - 1);
alpha_m = (2.5 - 0.1 * V_m) / (exp(2.5 -( 0.1 * V_m)) - 1);
alpha_h = (0.07) / (exp(0.05 * V_m));
beta_n = (0.125) / (exp( 0.0125 * V_m));
beta_m = (4)/(exp(V_m / 18));
beta_h = (1)/(exp(3-0.1 * V_m) + 1);

n(1,1) = alpha_n / (alpha_n + beta_n);
h(1,1) = alpha_h / (alpha_h + beta_h);
m(1,1) = alpha_m / (alpha_m + beta_m);
 
for i = 1:k-1
    
    alpha_n = (0.1 - (0.01 * V_m(i,1))) / (exp(1 - (0.1 * V_m(i,1))) - 1);
    alpha_m = (2.5 - 0.1 * V_m(i,1)) / (exp(2.5 -( 0.1 * V_m(i,1))) - 1);
    alpha_h = (0.07) / (exp(0.05 * V_m(i,1)));
    beta_n = (0.125) / (exp( 0.0125 * V_m(i,1)));
    beta_m = (4)/(exp(V_m(i,1) / 18));
    beta_h = (1)/(exp(3-(0.1 * V_m(i,1))) + 1);
    
    
    I_K = G_K_max * (n(i,1)^4) * (V_m(i,1) - V_K);
    I_Na = G_Na_max * (m(i,1)^3) * h(i,1) * (V_m(i,1) - V_Na);
    I_L = G_L_max * (V_m(i,1) - V_L);
    I_i = I_K + I_Na + I_L;
    I_m = I(i) - I_i;
    
    V_m(i,2) = V_m(i,1) + aralik * I_m / c_m;
    n(i,2) = n(i,1) + aralik * (alpha_n*(1-n(i,1))- beta_n*n(i,1));
    m(i,2) = m(i,1) + aralik * (alpha_m*(1-m(i,1))- beta_m*m(i,1));
    h(i,2) = h(i,1) + aralik * (alpha_h*(1-h(i,1))- beta_h*h(i,1));
    
    for j=2:p
                
        alpha_n = (0.1 - (0.01 * V_m(i,j))) / (exp(1 - (0.1 * V_m(i,j))) - 1);
        alpha_m = (2.5 - 0.1 * V_m(i,j)) / (exp(2.5 -( 0.1 * V_m(i,j))) - 1);
        alpha_h = (0.07) / (exp(0.05 * V_m(i,j)));
        beta_n = (0.125) / (exp( 0.0125 * V_m(i,j)));
        beta_m = (4)/(exp(V_m(i,j) / 18));
        beta_h = (1)/(exp(3-(0.1 * V_m(i,j))) + 1);
        
        I_Na = (m(i,j)^3) * G_Na_max * h(i,j) * (V_m(i,j-1) - V_Na);
        I_K = G_K_max * (n(i,j)^4) * (V_m(i,j-1) - V_K);
        I_L = G_L_max * (V_m(i,j-1) - V_L);
        I_i= I_K + I_Na + I_L;
        I_m = I(j) - I_i;
   
    
        V_m(i,j+1) = V_m(i,j) + aralik * I_m / c_m;
        n(i,j+1) = n(i,j) + aralik * (alpha_n*(1-n(i,j))- beta_n*n(i,j));
        m(i,j+1) = m(i,j) + aralik * (alpha_m*(1-m(i,j))- beta_m*m(i,j));
        h(i,j+1) = h(i,j) + aralik * (alpha_h*(1-h(i,j))- beta_h*h(i,j));
        
        
        v_m(i,j) = (aralik /(C_m*u_ar^2))*(V_m(i,j-1)-2*V_m(i,j)+V_m(i,j+1));
        V_m(i,j+1) = V_m(i,j) + v_m(i,j);
        
    end
    V_m(i+1,1) = V_m(i,1) + aralik * I_m / c_m;
    n(i+1,1) = n(i,1) + aralik * (alpha_n*(1-n(i,1))- beta_n*n(i,1));
    m(i+1,1) = m(i,1) + aralik * (alpha_m*(1-m(i,1))- beta_m*m(i,1));
    h(i+1,1) = h(i,1) + aralik * (alpha_h*(1-h(i,1))- beta_h*h(i,1));
    
    
end
    
V_m = V_m-90;

for r = 1:p
    V_m(1:k,r+1) = delayseq(V_m(1:k,r),(k-1)/10);
    V_m(1:(k-1)/10,r+1) = V_r ;
    
end
 
u = 1;
while u
    for n = 1:p
        plot(V_m(1:k,n),'LineWidth', 2)
        axis([0 k -140 40])
        title('Membran Voltajı')
        ylabel('Voltaj(mV)')
        xlabel('Akson Uzunlugu (cm)')
        pause(2)
    
        if n~=k
            clf
        end
        
    end
    u = u+1;
    if n == 11
        break
    end
    
    
end


%for j = 1:k   
%  plot(V_m(1:p,j))
%    title('Membran Potansiyeli')
%    ylabel('Voltaj (mv)')
%    xlabel('Akson uzunluğu (cm)')
%    pause(0.2);
%  
%end

%u = 1;

%while u
    
%    for j=1:k
%         plot(p(:,j),V_m);
%         hold on       
%         pause(0.01)
%         
%         if j~=k
%             clf
%         end
%    end
    
%    u=u+1;
%    if u==p
%        break
%    end
%end


%for j=1:5:p
%    plot(t,V_m(:,j),'b-');
%    title('Aksiyon Potansiyeli Yayılımı')
%    ylabel('Voltaj (mV)')
%    xlabel('Akson Uzunluğu (cm)')
%end


