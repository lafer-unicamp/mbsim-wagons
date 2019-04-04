files = list("freqtest-3d1Hz.mbsim.h5","freqtest-3d6Hz.mbsim.h5","freqtest-3d11Hz.mbsim.h5","freqtest-3d16Hz.mbsim.h5")

for f=files

// Open the created file
a = h5open(strcat(["~/git/mbsim-wagons/results/",f]));

d=h5read(a,"/Front truck/Bolster/C/data");
d = d';

plot(d(:,1),d(:,3),'color',rand(1,3));

// Free the resources
//h5close(a);
end;
