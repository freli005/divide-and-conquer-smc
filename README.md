# divide-and-conquer-smc
Code for the Ising model example in the paper 'Divide-and-Conquer with Sequential Monte Carlo'

This repository contains the matlab code used for the Ising model example (Section 5.1 of the main paper and Section C.1 of the supplementary material) for the paper:

  Divide-and-Conquer with Sequential Monte Carlo
  F. Lindsten, A. M. Johansen, C. A. Naesseth, B. Kirkpatrick, T. B. Schön, J. A. D. Aston, and A. Bouchard-Côté
  http://arxiv.org/abs/1406.4993

The parameters of the Ising model are specified in the function setup_model.m. The various samplers can be run using test_(name-of-sampler).m. Please bear in mind that this is research code and that we have experimented with various different settings for the D&C-SMC methodology. This is the reason for why a large number of parameters can be specified by the user, affecting the properties of the sampler in different ways as explained in the code.
