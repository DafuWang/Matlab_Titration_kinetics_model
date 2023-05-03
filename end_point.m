%%计算电导滴定的滴定液浓度、待测液体积最佳区域图，以及其它离子对最佳滴定区的影响
%Lambda代表溶液电导率，lambda是离子摩尔电导率
%Lambda_dt代表溶液电导率关于t的导数，lambda_dt是离子摩尔电导率关于t的导数

clear
clc
close all

ion=ion_lambda_G_z();%离子有关的电荷和稀溶液摩尔电导率及离子电导率经验系数

%滴定参数设定
upsilon=0.008e-3;%滴定速率
ion.SO4.C0_Na2SO4=0.01;%基准被滴定硫酸钠溶液设定
ion.Na.C0_Na2SO4=2*ion.SO4.C0_Na2SO4;%滴定液中的与硫酸根相应的Na
C0_BaCl2=0.005:0.005:0.1;%标准滴定溶液氯化钡浓度的设定，是一个范围也可以是一个值

%% 第1个循环是滴定液浓度的变化，第2个循环是其它离子的浓度的变化，第3个循环是待测液体积的变化
i=0;
for C_BaCl2=C0_BaCl2%
    i=i+1;
    ion.Ba.C0_BaCl2=C_BaCl2;%滴定浓度
    ion.Cl.C0_BaCl2=2*ion.Ba.C0_BaCl2;%滴定浓度

    %添加进去的影响离子浓度设定
    j=0;
    C0_bulk_NaCl=0.5:0.5:0.5;
    for ion_bulk_X=C0_bulk_NaCl
        j=j+1;
        ion.K.C0_bulk=0;
        ion.Ca.C0_bulk=0;
        ion.NO3.C0_bulk=0;

        ion.H.C0_bulk=0;
        ion.Cl.C0_bulk=ion_bulk_X;
        ion.Na.C0_bulk1=ion_bulk_X;
        ion.OH.C0_bulk=0;

        ion.Na.C0_bulk=ion.Na.C0_bulk1+ion.Na.C0_Na2SO4;%待测液中总的Na


        %离子初始参数总和计算，SBC分别代表S、Ba和Cl的影响
        [bulk,SBC]=ion_z_C0_G_lambda_infinity_bulk(ion);

        k=0;
        V_variation=0.1e-3:0.1e-3:30e-3;
        for V=V_variation%改变滴定溶液体积
            k=k+1;
            V3=V*ion.SO4.C0_Na2SO4/ion.Ba.C0_BaCl2;%理论滴定体积
            % %% 第一阶段
            % V1=0;%第一阶段结束的体积
            % t1=V1/upsilon;%达到第一阶段的时间
            % t_1=0:t1;%第一阶段的时间
            % [Lambda_stage1,Lambda_dt_stage1]=Stage1(bulk,SBC,V,upsilon,t_1);

            %% 第3阶段
            t3=V3/upsilon;%第3阶段结束的时间
            t_3=0:t3;%第3阶段的时间，忽略了第1和2阶段
            [Lambda_stage3,Lambda_dt_stage3,lambda_SO4_t0]=Stage3(bulk,SBC,V,upsilon,t_3);%lambda_SO4_t0产生的电导率
            Lambda_stage3_modified=Lambda_stage3;

            %% 第5阶段
            t_finty=1.2*V3/upsilon;%假设第5阶段结束的时间
            t_5=t3:t_finty;%第5阶段的时间，忽略了第4阶段

            [Lambda_stage5,Lambda_dt_stage5]=Stage5(bulk,SBC,V,upsilon,t3, t_5);
            Lambda_stage5_modified=Lambda_stage5+Lambda_stage3_modified(end)-Lambda_stage5(1);

            Volume_all=1000*upsilon*[t_3,t_5];
            Lambda_all=[Lambda_stage3_modified,Lambda_stage5_modified];
            X_Y{i,j,k}=[Volume_all',Lambda_all'];

            norm_Volume_all=Volume_all/(V3*1000);
            norm_Lambda_all=(Lambda_all-Lambda_stage3(end))/(abs(Lambda_stage3(1)-Lambda_stage3(end)));
            norm_X_Y{i,j,k}=[norm_Volume_all',norm_Lambda_all'];


            % figure(100*m+1);
            % plot(Volume_all,Lambda_all,'LineWidth',1.5)
            % set(0,'defaultfigurecolor','w');%设定背景颜色为白色
            % set(gca,'FontSize',10,'Fontname', 'Times New Roman');
            % xlabel('$Volume of titration (mL)$','Interpreter','latex','FontSize',15);
            % ylabel('$\Lambda (mS/cm)$','Interpreter','latex','FontSize',15);
            % set(gcf,'position',[360,198,560,420]);
            % set(gca,'position',[0.1,0.1,0.88,0.88]);
            % hold on

            % figure(n+3);
            % plot(Volume_all/(V3*1000),(Lambda_all-Lambda_stage3(end))/(abs(Lambda_stage3(end)-Lambda_stage3(1))),'LineWidth',1.5)
            % set(0,'defaultfigurecolor','w');%设定背景颜色为白色
            % set(gca,'FontSize',10,'Fontname', 'Times New Roman');
            % xlabel('$Volume of titration (mL)$','Interpreter','latex','FontSize',15);
            % ylabel('$\Lambda (mS/cm)$','Interpreter','latex','FontSize',15);
            % set(gcf,'position',[360,198,560,420]);
            % set(gca,'position',[0.1,0.1,0.88,0.88]);
            % hold on



            N=0;
            Volume_dt_all=1000*upsilon*[t_3(1:end-N),t_5(N+1:end)];
            Lambda_dt_all=[Lambda_dt_stage3(1:end-N),Lambda_dt_stage5(N+1:end)];
            X_Y_dt{i,j,k}=[Volume_dt_all',Lambda_dt_all'];

            norm_Volume_dt_all=Volume_dt_all/(V3*1000);
            norm_Lambda_dt_all=Lambda_dt_all*V3*1000/(abs(Lambda_stage3(end)-Lambda_stage3(1)));
            norm_X_Y_dt{i,j,k}=[norm_Volume_dt_all',norm_Lambda_dt_all'];

            
            % figure(10*m+2)
            % plot(Volume_dt_all,Lambda_dt_all,'LineWidth',1.5)
            % set(0,'defaultfigurecolor','w');%设定背景颜色为白色
            % set(gca,'FontSize',10,'Fontname', 'Times New Roman');
            % xlabel('$Volume of titration (mL)$','Interpreter','latex','FontSize',15);
            % ylabel('$\dot \Lambda(mS/(cm \cdot mL))$','Interpreter','latex','FontSize',15);
            % set(gcf,'position',[360,198,560,420]);
            % set(gca,'position',[0.11,0.1,0.87,0.86]);
            % hold on


            Lambda_dt_3(i,j,k)=Lambda_dt_stage3(end);
            Lambda_dt_5(i,j,k)=Lambda_dt_stage5(1);
            delta_Lambda_dt_5_3(i,j,k)=Lambda_dt_5(i,j,k)-Lambda_dt_3(i,j,k);

            norm_Lambda_dt_3(i,j,k)=Lambda_dt_stage3(end)*V3*1000/(abs(Lambda_stage3(end)-Lambda_stage3(1)));
            norm_Lambda_dt_5(i,j,k)=Lambda_dt_stage5(1)*V3*1000/(abs(Lambda_stage3(end)-Lambda_stage3(1)));
            norm_Lambda_dt_5_3(i,j,k)=norm_Lambda_dt_5(i,j,k)-norm_Lambda_dt_3(i,j,k);


          %  delta_SO4(i,j,k)=lambda_SO4_t0/Lambda_stage3(1,1);%初始溶液中硫酸根对待测液电导率的贡献

            theta_stage3(i,j,k)=atan(Lambda_dt_3(i,j,k))*180/pi;%未归一化滴定曲线滴定终点左切线夹角，角度值
            theta_stage5(i,j,k)=atan(Lambda_dt_5(i,j,k))*180/pi;%未归一化滴定曲线滴定终点右切线夹角，角度值

            norm_theta_stage3(i,j,k)=atan(norm_Lambda_dt_3(i,j,k))*180/pi;%归一化滴定曲线滴定终点左切线夹角，角度值
            norm_theta_stage5(i,j,k)=atan(norm_Lambda_dt_5(i,j,k))*180/pi;%归一化滴定曲线滴定终点右切线夹角，角度值


            delta_theta_5_3(i,j,k)=theta_stage5(i,j,k)-theta_stage3(i,j,k);%未归一化滴定曲线滴定终点左右切线夹角，角度值
            norm_delta_theta_5_3(i,j,k)=norm_theta_stage5(i,j,k)-norm_theta_stage3(i,j,k);%归一化滴定曲线滴定终点左右切线夹角，角度值

        end
    end
end


 
for i=1:size(C0_BaCl2,2)
    V_x=V_variation'*1000;%待测液体积的变化

    figure(7)%非归一化滴定曲线滴定终点左右切线的夹角随体积的变化
    y1(:,i)=delta_theta_5_3(i,j,:);
    plot(V_x,y1(:,i))
    hold on

    figure(8)%归一化滴定曲线滴定终点左右切线的夹角随体积的变化
    y2(:,i)=norm_delta_theta_5_3(i,j,:);
    plot(V_x,y2(:,i))
    hold on

end




% table=[Lambda_dt_3,Lambda_dt_5,delta_Lambda_dt_5_3,delta_theta_5_3,norm_Lambda_dt_5_3,norm_delta_theta_5_3];

% ylims=ylim;
% xlims=xlim;
% position=[0.8*xlims(1,2) 0.8*ylims(1,2) 0];
% set(text,'String',strcat('C(HCl)=',num2str(C_HCl),'mol/L'),'FontSize',10);















