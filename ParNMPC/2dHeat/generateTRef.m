% rng(2)
%% 
M = 13;
TRefRand = zeros(13,13);
idx = randi([1 13*13],1,1);
TRefRand(idx) = 1;
H = ones(9,9);
TRefAll = filter2(H,TRefRand)*100;
%% step
% TRefAll = zeros(13,13);
% TRefAll(1:6,:) = 200;
%% gradient
TRefAll = zeros(13,13);
for i=1:13
    TRefAll(:,i) = i/12*300-25;
end
%% V1
% TRefAll = zeros(13,13);
% for i=8:13
%     TRefAll(i,:) =  (i-8)/5*300;
% end
% for i=1:6
%     TRefAll(i,:) =  300-(i)/5*300+60;
% end
%% V2
TRefAll = zeros(13,13);
for i=9:13
    TRefAll(i,:) =  (i-9)/4*300;
end
for i=1:5
    TRefAll(i,:) =  300-(i)/4*300 + 75;
end
%% block
% TRefAll = zeros(13,13);
% TRefAll(5:9,5:9) = 400;


%% plot
surf(TRefAll+300);%heatmap
colormap jet
zlim([300 700]);
caxis([300 700]);
% view(0,90)
shading interp

drawnow


save TRefAll.mat TRefAll