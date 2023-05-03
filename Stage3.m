function [Lambda_stage3,Lambda_dt_stage3,lambda_SO4_t0]=Stage3(bulk,SBC,V,upsilon,t_3)
% 计算除了硫酸根和滴定离子外其他离子对离子强度的贡献,及其导数,每个时间对应一个数
N=size(bulk.z,2);
    Gamma_stage3=0;Gamma_dt_stage3=0;
for i=1:N
    Gamma_stage3=   Gamma_stage3    + 0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_3 );
    Gamma_dt_stage3=Gamma_dt_stage3 + 0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_3 ).^2.*(-upsilon);
end
%第2项
    Gamma_stage3=   Gamma_stage3    + 0.5 *  SBC.z(1)^2 .* V * SBC.C0(1) ./ (V + upsilon .* t_3);
    Gamma_dt_stage3=Gamma_dt_stage3 + 0.5 *  SBC.z(1)^2 .* V * SBC.C0(1) ./ (V + upsilon .* t_3).^2.*(-upsilon);

%第3项 
    Gamma_stage3=   Gamma_stage3    - 0.5 *  SBC.z(1)^2 .* upsilon .* SBC.C0(2) .* t_3 ./ (V + upsilon .* t_3);
    Gamma_dt_stage3=Gamma_dt_stage3 - 0.5 *  SBC.z(1)^2 .* upsilon .* SBC.C0(2) .* V ./ (V + upsilon .* t_3).^2;
%第4项
    Gamma_stage3=   Gamma_stage3    + 0.5 *  SBC.z(3)^2 .* upsilon .* SBC.C0(3) .* t_3 ./ (V + upsilon .* t_3);
    Gamma_dt_stage3=Gamma_dt_stage3 + 0.5 *  SBC.z(3)^2 .* upsilon .* SBC.C0(3) .* V ./ (V + upsilon .* t_3).^2;

%% 计算摩尔电导率及其导数，每种离子每个时间对应一个数
for i=1:N
lambda_bulk_stage3(i,:)=bulk.lambda_infinity(i) ./ (1+ bulk.G(i) .* sqrt(Gamma_stage3));
lambda_bulk_dt_stage3(i,:)=-0.5*bulk.lambda_infinity(i).*bulk.G(i).*Gamma_dt_stage3./(sqrt(Gamma_stage3).*(1+bulk.G(i).*sqrt(Gamma_stage3)).^2);
end
%第2项
i=i+1;j=1;
lambda_SBC_stage3(i,:)= SBC.lambda_infinity(j) ./ (1+  SBC.G(j) .* sqrt(Gamma_stage3));
lambda_SBC_dt_stage3(i,:)=-0.5*SBC.lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage3./(sqrt(Gamma_stage3).*(1+SBC.G(j).*sqrt(Gamma_stage3)).^2);
lambda_SO4_t0=lambda_SBC_stage3(i,1);
%第3项 第3阶段第3项用到的稀溶液参数仍然是硫酸根的
i=i+1;j=2;
lambda_SBC_stage3(i,:)=SBC.lambda_infinity(j-1) ./ (1+ SBC.G(j-1) .* sqrt(Gamma_stage3));
lambda_SBC_dt_stage3(i,:)=-0.5*SBC.lambda_infinity(j-1).*SBC.G(j-1).*Gamma_dt_stage3./(sqrt(Gamma_stage3).*(1+SBC.G(j-1).*sqrt(Gamma_stage3)).^2);
%第4项
i=i+1;j=3;
lambda_SBC_stage3(i,:)=SBC.lambda_infinity(j) ./ (1+ SBC.G(j) .* sqrt(Gamma_stage3));
lambda_SBC_dt_stage3(i,:)=-0.5*SBC. lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage3./(sqrt(Gamma_stage3).*(1+SBC.G(j).*sqrt(Gamma_stage3)).^2);
%% 计算溶液电导率及其导数，每个时间对应一个数
% 计算除了硫酸根和滴定离子外所设计的其他参数，
a=0;
Lambda_stage3=0;Lambda_dt_stage3=0;
for i=1:N
    Lambda_stage3=Lambda_stage3 + bulk.z(i) .* lambda_bulk_stage3(i,:) .* V * bulk.C0(i) ./ (V + upsilon .* t_3);
   
    a0=bulk.z(i) .* V * bulk.C0(i);  
    b0=(V + upsilon .* t_3).^2; 
    c0=lambda_bulk_dt_stage3(i,:).*(V + upsilon .* t_3)-lambda_bulk_stage3(i,:).*upsilon;
    Lambda_dt_stage3=Lambda_dt_stage3+a0.*c0./b0;
end
%第2项参数
i=i+1;j=1;
Lambda_stage3=Lambda_stage3 + SBC.z(j) .* lambda_SBC_stage3(i,:) .* V * SBC.C0(j) ./ (V + upsilon .* t_3);

a1=SBC.z(j) .* V * SBC.C0(j);  
b1=(V + upsilon .* t_3).^2; 
c1=lambda_SBC_dt_stage3(i,:).*(V + upsilon .* t_3)-lambda_SBC_stage3(i,:).*upsilon;
Lambda_dt_stage3=Lambda_dt_stage3+a1.*c1./b1;
%第3项参数,阶段3的第3项仅离子浓度用的是Ba的，其他离子z，SBC.C0(2)，lambda率和G用的是SO4，%相对于第一阶段更改了成了减号，等式前是减号
i=i+1;j=2;
Lambda_stage3=Lambda_stage3 - SBC.z(j-1) .* lambda_SBC_stage3(i,:) .* upsilon .* t_3 * SBC.C0(j) ./ (V + upsilon .* t_3);

a2=SBC.z(j-1) .* upsilon.* SBC.C0(j);
b2=lambda_SBC_dt_stage3(i,:).*t_3+lambda_SBC_stage3(i,:);
c2=V + upsilon .* t_3;
d2=lambda_SBC_stage3(i,:).*upsilon.*t_3;
e2=(V + upsilon .* t_3).^2; 
Lambda_dt_stage3=Lambda_dt_stage3 - a2.*(b2.*c2-d2)./e2;
%第4项参数
i=i+1;j=3;
Lambda_stage3=Lambda_stage3 + SBC.z(j) .* lambda_SBC_stage3(i,:) .* upsilon .* t_3 * SBC.C0(j) ./ (V + upsilon .* t_3);

a3=SBC.z(j) .* upsilon.* SBC.C0(j);
b3=lambda_SBC_dt_stage3(i,:).*t_3 + lambda_SBC_stage3(i,:);
c3=V + upsilon .* t_3;
d3=lambda_SBC_stage3(i,:).*upsilon.*t_3;
e3=(V + upsilon .* t_3).^2; 
Lambda_dt_stage3=Lambda_dt_stage3+a3.*(b3.*c3-d3)./e3;
Lambda_dt_stage3=Lambda_dt_stage3./upsilon/1000;








