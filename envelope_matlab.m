% Test Swanepoel

num = xlsread('espectro.xlsx');
lambdas = num(250:end,1);
film = num(250:end, 2);
scatter(lambdas, film)
[up, lo] = envelope(film, 500);
hold on
plot(lambdas, up, lambdas, lo)
hold off