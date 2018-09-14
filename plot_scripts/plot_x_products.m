std_Y = 2.1505;
t = linspace(-3,3,20)';
Y = zeros(length(t));
for ii = 1:length(t)
    for jj = 1:length(t)
    	Y(ii,jj) =  generate_Y_values([t(ii),t(jj)],2,1,0,std_Y);
    end
end

surf(t,t,Y)
grid on
colormap jet
colorbar