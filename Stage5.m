function [Lambda_stage5,Lambda_dt_stage5]=Stage5(bulk,SBC,V,upsilon,t3, t_5)
% 计算除了硫酸根和滴定离子外其他离子对离子强度的贡献,及其导数,每个时间对应一个数
N=size(bulk.z,2);
    Gamma_stage5=0;Gamma_dt_stage5=0;
for i=1:N
    Gamma_stage5=   Gamma_stage5    + 0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_5 );
    Gamma_dt_stage5=Gamma_dt_stage5 + 0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_5 ).^2.*(-upsilon);
end
%第2项 upsilon.*t3是常数，相当于之前的V
    j=2;
    Gamma_stage5=   Gamma_stage5    - 0.5 *  SBC.z(j)^2 .* upsilon.*t3 * SBC.C0(j) ./ (V + upsilon .* t_5);
    Gamma_dt_stage5=Gamma_dt_stage5 - 0.5 *  SBC.z(j)^2 .* upsilon.*t3 * SBC.C0(j) ./ (V + upsilon .* t_5).^2.*(-upsilon);

%第3项 
    j=2;
    Gamma_stage5=   Gamma_stage5    + 0.5 *  SBC.z(j)^2 .* upsilon .* SBC.C0(j) .* t_5 ./ (V + upsilon .* t_5);
    Gamma_dt_stage5=Gamma_dt_stage5 + 0.5 *  SBC.z(j)^2 .* upsilon .* SBC.C0(j) .* V ./ (V + upsilon .* t_5).^2;
%第4项
    j=3;
    Gamma_stage5=   Gamma_stage5    + 0.5 *  SBC.z(j)^2 .* upsilon .* SBC.C0(j) .* t_5 ./ (V + upsilon .* t_5);
    Gamma_dt_stage5=Gamma_dt_stage5 + 0.5 *  SBC.z(j)^2 .* upsilon .* SBC.C0(j) .* V ./ (V + upsilon .* t_5).^2;

%% 计算摩尔电导率及其导数，每种离子每个时间对应一个数
%第5阶段就和SO4无关了
for i=1:N
lambda_bulk_stage5(i,:)=bulk.lambda_infinity(i) ./ (1+ bulk.G(i) .* sqrt(Gamma_stage5));
lambda_bulk_dt_stage5(i,:)=-0.5*bulk.lambda_infinity(i).*bulk.G(i).*Gamma_dt_stage5./(sqrt(Gamma_stage5).*(1+bulk.G(i).*sqrt(Gamma_stage5)).^2);
end
%第2项 第5阶段第2项用到的都是Ba了
i=i+1;j=2;
lambda_SBC_stage5(i,:)= SBC.lambda_infinity(j) ./ (1+  SBC.G(j) .* sqrt(Gamma_stage5));
lambda_SBC_dt_stage5(i,:)=-0.5*SBC.lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage5./(sqrt(Gamma_stage5).*(1+SBC.G(j).*sqrt(Gamma_stage5)).^2);

%第3项 第5阶段第2项用到的也都是Ba
i=i+1;j=2;
lambda_SBC_stage5(i,:)=SBC.lambda_infinity(j) ./ (1+ SBC.G(j) .* sqrt(Gamma_stage5));
lambda_SBC_dt_stage5(i,:)=-0.5*SBC.lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage5./(sqrt(Gamma_stage5).*(1+SBC.G(j).*sqrt(Gamma_stage5)).^2);
%第4项
i=i+1;j=3;
lambda_SBC_stage5(i,:)=SBC.lambda_infinity(j) ./ (1+ SBC.G(j) .* sqrt(Gamma_stage5));
lambda_SBC_dt_stage5(i,:)=-0.5*SBC. lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage5./(sqrt(Gamma_stage5).*(1+SBC.G(j).*sqrt(Gamma_stage5)).^2);
%% 计算溶液电导率及其导数，每个时间对应一个数
% 计算除了硫酸根和滴定离子外所设计的其他参数，
a=0;
Lambda_stage5=0;Lambda_dt_stage5=0;
for i=1:N
    Lambda_stage5=Lambda_stage5 + bulk.z(i) .* lambda_bulk_stage5(i,:) .* V * bulk.C0(i) ./ (V + upsilon .* t_5);
   
    a0=bulk.z(i) .* V * bulk.C0(i);  
    b0=(V + upsilon .* t_5).^2; 
    c0=lambda_bulk_dt_stage5(i,:).*(V + upsilon .* t_5)-lambda_bulk_stage5(i,:).*upsilon;
    Lambda_dt_stage5=Lambda_dt_stage5+a0.*c0./b0;
end
%第2项参数 阶段5的第2项只与Ba有关，%相对于第一阶段更改了成了减号，等式前是减号
i=i+1;j=2;
Lambda_stage5=Lambda_stage5 - SBC.z(j) .* lambda_SBC_stage5(i,:) .* upsilon.*t3 .* SBC.C0(j) ./ (V + upsilon .* t_5);

a1=SBC.z(j) .* upsilon.*t3 * SBC.C0(j);  
b1=(V + upsilon .* t_5).^2; 
c1=lambda_SBC_dt_stage5(i,:).*(V + upsilon .* t_5) - lambda_SBC_stage5(i,:).*upsilon;
Lambda_dt_stage5=Lambda_dt_stage5 - a1.*c1./b1;
%第3项参数, 阶段5的第3项只与Ba有关
i=i+1;j=2;
Lambda_stage5=Lambda_stage5 + SBC.z(j) .* lambda_SBC_stage5(i,:) .* upsilon .* t_5 * SBC.C0(j) ./ (V + upsilon .* t_5);

a2=SBC.z(j) .* upsilon.* SBC.C0(j);
b2=lambda_SBC_dt_stage5(i,:).*t_5+lambda_SBC_stage5(i,:);
c2=V + upsilon .* t_5;
d2=lambda_SBC_stage5(i,:).*upsilon.*t_5;
e2=(V + upsilon .* t_5).^2; 
Lambda_dt_stage5=Lambda_dt_stage5 + a2.*(b2.*c2-d2)./e2;
%第4项参数
i=i+1;j=3;
Lambda_stage5=Lambda_stage5 + SBC.z(j) .* lambda_SBC_stage5(i,:) .* upsilon .* t_5 * SBC.C0(j) ./ (V + upsilon .* t_5);

a3=SBC.z(j) .* upsilon.* SBC.C0(j);
b3=lambda_SBC_dt_stage5(i,:).*t_5 + lambda_SBC_stage5(i,:);
c3=V + upsilon .* t_5;
d3=lambda_SBC_stage5(i,:).*upsilon.*t_5;
e3=(V + upsilon .* t_5).^2; 
Lambda_dt_stage5=Lambda_dt_stage5+a3.*(b3.*c3-d3)./e3;
Lambda_dt_stage5=Lambda_dt_stage5./upsilon/1000;








