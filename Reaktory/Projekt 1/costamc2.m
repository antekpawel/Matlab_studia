function costamc2(x,jeden,dwa,c)
figure(x)
clf
hold all
grid on
grid minor
ylabel('$$ {Q \over Q_{max}} $$','Interpreter','latex')
xlabel('$$ {L \over L_{max}} $$','Interpreter','latex')
title(' ')
plot(jeden(:,3),jeden(:,c),'r-x')
plot(dwa(:,3),dwa(:,c),'b-x')
legend('Otwieranie','Zamykanie','Location','Northwest')
end