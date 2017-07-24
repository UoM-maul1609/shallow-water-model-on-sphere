% do the fourier analysis
%fileNames={'output_0_2.nc','output_1_2.nc','output_2_2.nc','output_3_2.nc',...
%    'output_4_2.nc',...
%    'output_5_2.nc','output_6_2.nc','output_7_2.nc','output_8_2.nc',...
%    'output_9_2.nc','output_10_2.nc','output_11_2.nc','output_12_2.nc',...
%    'output_13_2.nc','output_14_2.nc','output_15_2.nc','output_16_2.nc',...
%    'output_17_3.nc'};

%u_jet=[5. 10. 20. 30. 40. 50. 60. 70. 80. 90. 100. 125. 150. 175. 200. 250 300 350];

fileNames={'/tmp/output.nc'};
u_jet=[100.];
figure;
fourier_wave_number(fileNames,u_jet);
normal_modes_compare(u_jet);
title('Model (circles) vs Linear theory (lines)');