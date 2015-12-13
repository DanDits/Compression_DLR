# Compression_DLR
Compression using a dynamic low rank scheme for matrices and matrix differential equations

Based on paper "A projector-splitting integrator for dynamical low-rank approximation" from Christian Lubich and Ivan V. Oseledets published in Springer in 2013.

Implements the general DLR Routine with the core being the DLR_Step. 
This routine is used 
1. In DLR_Frames on a list of Frames that come from some generated movies or gif files.
  2. In DLR_Integrator to demonstrate the integration capability of the method by applying it to the inhomogenous heat equation with some arbitrary time and position dependend inhomogenity f.
  
The compression level is controlled by the rank parameter. This can either be a value between 0.0 and 2.0 with 0. being maximal compression (rank=1), 1.0 being no compression (rank = max dimension * 0.5) and 2.0 being under compression (rank = max dimension). The rank can also be supplied directly as a number. The rank approximation can be chosen to be adaptively or fixed. 

This does NOT produce any official file format for the movies but only offers the matrices U,S,V with U*S*V'=Frame. The compression happens since U and V both have only r columns and S is an rxr matrix with r being the chosen compression rank. (See singular value decomposition for the concept). This offers a fast calculation of these decompositions of a series of frames or a time dependend differential equation of order 1.
