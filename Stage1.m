function [Lambda_stage1,Lambda_dt_stage1]=Stage1(bulk,SBC,V,upsilon,t_1)

% 计算除了硫酸根和滴定离子外其他离子对离子强度的贡献,及其导数,每个时间对应一个数
N=size(bulk.z,2);
    Gamma_stage1=0;Gamma_dt_stage1=0;
for i=1:N
    Gamma_stage1=   Gamma_stage1 +    0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_1 );
    Gamma_dt_stage1=Gamma_dt_stage1 + 0.5 * bulk.z(i)^2. * V * bulk.C0(i) ./ (V + upsilon .*t_1 ).^2.*(-upsilon);
end
%第2项
    Gamma_stage1=   Gamma_stage1 +    0.5 *  SBC.z(1)^2 .* V * SBC.C0(1) ./ (V + upsilon .* t_1);
    Gamma_dt_stage1=Gamma_dt_stage1 + 0.5 *  SBC.z(1)^2 .* V * SBC.C0(1) ./ (V + upsilon .* t_1).^2.*(-upsilon);

%第3项
    Gamma_stage1=   Gamma_stage1 +    0.5 *  SBC.z(2)^2 .* upsilon .* SBC.C0(2) .* t_1 ./ (V + upsilon .* t_1);
    Gamma_dt_stage1=Gamma_dt_stage1 + 0.5 *  SBC.z(2)^2 .* upsilon .* SBC.C0(2) .* V ./ (V + upsilon .* t_1).^2;
%第4项
    Gamma_stage1=   Gamma_stage1 +    0.5 *  SBC.z(3)^2 .* upsilon .* SBC.C0(3) .* t_1 ./ (V + upsilon .* t_1);
    Gamma_dt_stage1=Gamma_dt_stage1 + 0.5 *  SBC.z(3)^2 .* upsilon .* SBC.C0(3) .* V ./ (V + upsilon .* t_1).^2;

%% 计算摩尔电导率及其导数，每种离子每个时间对应一个数
for i=1:N
lambda_bulk_stage1(i,:)=bulk.lambda_infinity(i) ./ (1+ bulk.G(i) .* sqrt(Gamma_stage1));
lambda_bulk_dt_stage1(i,:)=-0.5*bulk.lambda_infinity(i).*bulk.G(i).*Gamma_dt_stage1./(sqrt(Gamma_stage1).*(1+bulk.G(i).*sqrt(Gamma_stage1)).^2);
end
%第2项
i=i+1;j=1;
lambda_SBC_stage1(i,:)= SBC.lambda_infinity(j) ./ (1+  SBC.G(j) .* sqrt(Gamma_stage1));
lambda_SBC_dt_stage1(i,:)=-0.5*SBC.lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage1./(sqrt(Gamma_stage1).*(1+SBC.G(j).*sqrt(Gamma_stage1)).^2);

%第3项
i=i+1;
lambda_SBC_stage1(i,:)=SBC.lambda_infinity(j) ./ (1+ SBC.G(j) .* sqrt(Gamma_stage1));
lambda_SBC_dt_stage1(i,:)=-0.5*SBC.lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage1./(sqrt(Gamma_stage1).*(1+SBC.G(j).*sqrt(Gamma_stage1)).^2);
%第4项
i=i+1;
lambda_SBC_stage1(i,:)=SBC.lambda_infinity(j) ./ (1+ SBC.G(j) .* sqrt(Gamma_stage1));
lambda_SBC_dt_stage1(i,:)=-0.5*SBC. lambda_infinity(j).*SBC.G(j).*Gamma_dt_stage1./(sqrt(Gamma_stage1).*(1+SBC.G(j).*sqrt(Gamma_stage1)).^2);
%% 计算溶液电导率及其导数，每个时间对应一个数
% 计算除了硫酸根和滴定离子外所设计的其他参数，
a=0;
Lambda_stage1=0;Lambda_dt_stage1=0;
for i=1:N
    Lambda_stage1=Lambda_stage1 + bulk.z(i) .* lambda_bulk_stage1(i,:) .* V * bulk.C0(i) ./ (V + upsilon .* t_1);
   
    a0=bulk.z(i) .* V * bulk.C0(i);  
    b0=(V + upsilon .* t_1).^2; 
    c0=lambda_bulk_dt_stage1(i,:).*(V + upsilon .* t_1)-lambda_bulk_stage1(i,:).*upsilon;
    Lambda_dt_stage1=Lambda_dt_stage1+a0.*c0./b0;
end
%第2项参数
i=i+1;j=1;
Lambda_stage1=Lambda_stage1 + SBC.z(j) .* lambda_SBC_stage1(i,:) .* V * SBC.C0(j) ./ (V + upsilon .* t_1);

a1=SBC.z(j) .* V * SBC.C0(j);  
b1=(V + upsilon .* t_1).^2; 
c1=lambda_SBC_dt_stage1(i,:).*(V + upsilon .* t_1)-lambda_SBC_stage1(i,:).*upsilon;
Lambda_dt_stage1=Lambda_dt_stage1+a1.*c1./b1;
%第3项参数
i=i+1;j=2;
Lambda_stage1=Lambda_stage1 + SBC.z(j) .* lambda_SBC_stage1(i,:) .* upsilon .* t_1 * SBC.C0(j) ./ (V + upsilon .* t_1);

a2=SBC.z(j) .* upsilon.* SBC.C0(j);
b2=lambda_SBC_dt_stage1(i,:).*t_1+lambda_SBC_stage1(i,:);
c2=V + upsilon .* t_1;
d2=lambda_SBC_stage1(i,:).*upsilon.*t_1;
e2=(V + upsilon .* t_1).^2; 
Lambda_dt_stage1=Lambda_dt_stage1+a2.*(b2.*c2-d2)./e2;
%第4项参数
i=i+1;j=3;
Lambda_stage1=Lambda_stage1 + SBC.z(j) .* lambda_SBC_stage1(i,:) .* upsilon .* t_1 * SBC.C0(j) ./ (V + upsilon .* t_1);

a3=SBC.z(j) .* upsilon.* SBC.C0(j);
b3=lambda_SBC_dt_stage1(i,:).*t_1 + lambda_SBC_stage1(i,:);
c3=V + upsilon .* t_1;
d3=lambda_SBC_stage1(i,:).*upsilon.*t_1;
e3=(V + upsilon .* t_1).^2; 
Lambda_dt_stage1=Lambda_dt_stage1+a3.*(b3.*c3-d3)./e3;
Lambda_dt_stage1=Lambda_dt_stage1./upsilon/1000;







