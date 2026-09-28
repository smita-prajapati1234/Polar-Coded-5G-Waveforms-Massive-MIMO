clear all
clc
%%------------------------------------------------------------------------------------------------------------%%
%                                           Design of Prototype Filter 
%%------------------------------------------------------------------------------------------------------------%%
M=16;% number of sub-channels
%% Prototype Filter For FBMC 
K=4; % overlapping factor 
lp=K*M-1;% prototype filter length
y=[1 0.97195983 sqrt(2)/2 0.23514695 ]; % coefficient of filter when K=4
yy=[1 0.911438 0.411438];% coefficient of filter when K=3
yyy=[1 sqrt(2)/2];% coefficient of filter when K=2
u=lp-1;
s=2*pi/(K*M);
for m=0:u-1
    r=m+1;
p1(m+1)=y(1,1)-2*y(1,2)*cos(r*s)+2*y(1,3)*cos(2*r*s)-2*y(1,4)*cos(3*r*s);% prototype filter equation K=4
p2(m+1)=yy(1,1)-2*yy(1,2)*cos(r*s)+2*yy(1,3)*cos(2*r*s); % prototype filter equation K=3
p3(m+1)=yyy(1,1)-2*yyy(1,2)*cos(r*s); % prototype filter equation K=2
end
%% Prototype Filter For OFDM
for m=0:M-1
p4(m+1)=1;
end
%% Frequency response of prototype filter for FBMC and OFDM
fvtool(p1,1,p4,1)
legend('FBMC','OFDM')
fvtool(p1,1,p2,1,p3,1)
legend('K=4','K=3','K=2')
