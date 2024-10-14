% testing data from
% https://www.sciencedirect.com/science/article/pii/S0896844615000716?via%3Dihubs
data=readmatrix('data2.csv');
exp23=data(data(:,2)==23,:);
exp33=data(data(:,2)==33,:);
exp33=exp33(exp33(:,1)>=50,:);

[E23,S23,I23]=activation_energy(exp23(:,1)+272.15,exp23(:,3));
[E33,S33,I33]=activation_energy(exp33(:,1)+272.15,exp33(:,3));

figure, hold on
plot(1./(exp23(:,1)+272.15),exp23(:,3),'k.')
plot(1./(exp33(:,1)+272.15),exp33(:,3),'r.')
plot(1./(exp23(:,1)+272.15),I23+S23./(exp23(:,1)+272.15),'k')
plot(1./(exp33(:,1)+272.15),I33+S33./(exp33(:,1)+272.15),'r')

function [Ea,slope,intercept]=activation_energy(temp,logK)
% Rearranges Arrhenius Equation to caclulate activation energy (kJ/mol)
% from temperature(K) and log(dissolution rate)

inv_temp=1./temp;
X=[ones(length(inv_temp),1) inv_temp];
mc=X\logK;

Ea=-(8.314*log(10))*mc(2)*1e-3; % Activation energy in kJ/mol (Log(10) to convert from ln to log)
intercept=mc(1);
slope=mc(2);
end
