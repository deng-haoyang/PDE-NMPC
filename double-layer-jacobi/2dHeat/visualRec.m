M    = 13;
TMIN = 300;
TMAX = 700;

%%
data = load('GEN_log_rec.mat');
rec = data.rec;
T = zeros(M,M);
recSize = length(rec.x(:,1));

for tStep = 1:1:recSize-1%
    uCnt = 1;
    xCnt = 1;
    for yi=1:M
        for xi=1:M
            if mod(xi-1,4)==0 && mod(yi-1,4)==0
                T(xi,yi) = rec.u(tStep,uCnt);
                uCnt     = uCnt + 1;
            else
                T(xi,yi) = rec.x(tStep+1,xCnt);
                xCnt     = xCnt + 1;
            end
        end
    end
    
%   
    figure(2);
%     surf(T,'EdgeColor','interp','FaceColor','interp');%heatmap
    surf(T);
%     hold on
    colormap jet	

    shading interp

    zlim([TMIN TMAX]);
    caxis([TMIN TMAX]);
%      view(0,90)

%     
%     heatmap(T);
%     colormap jet
%     caxis([300 350]);
    
    
    drawnow
    tStep
    pause(0.1);

    
end
