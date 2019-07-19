%% Frequency and Coupling Constant 개형 보기
mode_n=5;

spectral_function = figure;
plot(W/eV,-2*imag(G))
xlabel('\omega/eV')
ylabel('A(\omega)')
title(sprintf('250000 sqrtN-1 qTF = 2, mode = %g, broadening 1e-5 ',mode_n ))
saveas(spectral_function,sprintf('250000 sqrtN-1 qTF = 2,  mode = %g, broadening 1e-5.svg',mode_n))