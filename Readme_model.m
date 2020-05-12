%% To simulate models for resting state activity using
%% Spiking and DMF mean field model at each cortical node
%% along with Synaptic Plasticity model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Resting state activity simulation using DMF model
%%  and generating BOLD activity
%% To generate figure 8 in the new manuscript need to
%% obtain bifurcation diagram firing rate use 
%% Eq.(7) mean field model section in the manuscript

DMF_model_bold_balance
%% you also need to estmate excitatory firing rate conversion using sigmoidal function
phie
phii
%% To calculate BOLD signal from simulated neuronal activity
%% use 
BOLD
%% For BOLD plotting 
%% use 
Boldcalcandplot
%% For generating a figure like the one in manuscript 
%% you need small a SC based on [3-6] cortical nodes that
%% are implemented
%% Provided here is SC_capacity_normalized.mat capacity as well as distance matrix
%% and also simulated FC_rest, FC_plastic

plotfigure_rsFCcorrelations.m

%% for a custom made colormap we have used in our article
%% use the following script
darkb2r

%% To simulate full spiking network mean field model in the 
%% absence of plasticity use a C environment as it is computationally
%% expensive  
%% main manuscript these equations appear in Eq.(1-5)

%% based on Deco simulation of mean field activity
%% to obtain bifurcatin digram as a function of global coupling
%% Cortical pool is scripted in 

%%%pool.cpp
%%%pool.h

%% Pool related computations are scripted in

%%%Brain.cpp
%%%Brain.h

%% finally mean field simulation can be carried out using 
%% Global connectivity SC matrix
%%%mfgus.cpp

%% STDP is implemented in a standalone Python environment
%% Brian spiking network is used to construct pool activity
%% Subsequently relevant node model neurons are connected 
%% short range also long range synapses.
%% It simulates full network with brain state dependent
%% input and exponential STDP
%% For simplification from node level synaptic dynamics 
%% numerically estimated using both in the presence as well
%% as absence of NMDA synapses.
%% Without NMDA synapses
%% Main manuscript these Eqs appear in (1-5) and also in Eq.(8-9)

%%%DRstdpnode_model2013_BrianCode.py

%% With NMDA synapses 
%%%DR_stdpmodel_withNMDA.py

%% output spiking activity, Network conductances, Mean population activity 
%% Modification of local excitatory-excitatory weights.
%% is monitored and plotted in a separate script

%%%plot_paperfig4and5.py

%% For further analysis mean activity of network is saved 
%% all network activity are saved in MAT files.  
%% Time frequency analysis is carried out using multiple 
%% scripts

%%%Wavelet_1ch.m %% this script also used to plot for single subject single channel EEG

%%%%% main script that calculates wavelet transform
%%%wavelet33.m
%% for calulating fourier basis vectors 
%%%wave_bases.m 
%% finally for plotting a figure like the one 
%% manuscript labeled as figure 6
%% First load network pooled activity 
%% This script first calculates amplitude based on the
%% Low to high frequency resolution in the monitored network activity

%%%figure_plot_power_frequency.m 




