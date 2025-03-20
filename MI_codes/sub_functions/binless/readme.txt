Binned method notes:

- BlahutArimoto.m is used to find channel capacity if the distribution of signals is not specified in bin_CC_MI_js.m
- equaldensitybins.m is used to get bins with equal densities, which is supposed to reduce error in calculation of MI and CC


Binless method notes:

 - getCC.m, which taes an input cell vector. Each cell corresponds to different experimental conditions (such as different concentrations). Within each cell you have MxN matrix, where:
 	- M corresponds to number of measurements (1 if you are looking at max response, or more if you use multiple measurements or whole timeseries)
 	- N is the number of samples that you have. 
 The output is info in bits, and the optimal distribution of the input cell vector. In the inputs, the second input value is 1 or 0 depending on whether you want the fmincon optimization readout. The third possible input is if you want to specify the cell vector distribution,s in case you want to measure mutual information rather than channel capacity

- jacknifeCC extrapolates CC at an infinite number of data points (running getCC several times for different subsets of the data).
The input is the same, and the output is [fitCC,fullCC,I,Q], where fitCC is the jacknife calculation, 
fullCC is the standard getCC output, I and Q are all CC calculations that it used to get fitCC



which screws up the density estimation based on nearest neighbour distances. There is a requirement for the distribution of responses to be continuous, which in your case it is not around that zero value. I had that line in jacknifeCC so that you don't have repeated values, by adding some small noise to the data. In your data, it appears that the amount of noise needed is more than just eps that I specified in jacknife. By adding that noise, we are lowering the amount of information however, so there might be couple of ways to go about it. One is to remove all zero values from your calculation. The other is to balance adding noise to circumvent the continuous distribution requirement with lowering of the channel capacity

as the amount of added noise is reduced down to 1, we see a significant increase in CC. Subsequent reduction appears to incrementally increase CC, most likely due to inability to estimate density properly.