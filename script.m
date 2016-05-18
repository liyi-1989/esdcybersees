mylon=[-116.25,-113.75];
mylat=[43,45];
[Y,A]=get_obs(mylon,mylat);
csvwrite('obs_avg.csv',Y);