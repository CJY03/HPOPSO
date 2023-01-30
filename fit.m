function [ f ] = fit(x)
global h Wssum Wcsum Vck Eck Time;
%%%%%%动态线损
i=1;
 e(i)=(Wssum(i)-Wcsum(i)-x(2)-x(3)-Eck(i)*Vck*x(1))^2;
for i=2:h
    e(i)=(Wssum(i)-Wcsum(i)-x(2)-x(3)*((Wssum(i)^2*Time(1))/(Wssum(1)^2*Time(i)))-Eck(i)*Vck*x(1))^2;
end
%%%%%%固定线损
% for i=1:h
%  e(i)=(Wssum(i)-Wcsum(i)-x(2)-x(3)-Eck(i)*Vck*x(1))^2;
% end

%%%%%%无线损
% for i=1:h
%  e(i)=(Wssum(i)-Wcsum(i)-Eck(i)*Vck*x(1))^2;
% end

f=0;
for i=1:h
    f=f+e(i);
end

end


