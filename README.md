# code for 'Quickest Causal Change Point Detection by Adaptive Intervention'

main.m is the main file, which provides the parameter settings for the different synthetic experiments and two real case experiments discussed in the paper ‘Quickest Causal Change Point Detection by Adaptive Intervention’. After selecting the appropriate parameters, running main.m will generate the corresponding EDD v.s. ARL plot for the MAX-AI, MULTI-AI methods along with the other four baselines.

getW.m corresponds to our algorithm, incorporating different monitoring and intervention strategies. getW_maxOra.m and getW_multiOra.m correspond to the two oracle methods, respectively, while getW_GGM.m and getW_JSWL.m correspond to the other two baseline methods.
