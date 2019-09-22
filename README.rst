.. |date| date::

**********************************************
PPTC PDX Classifier and Correlation Figures
**********************************************

:author: Jo Lynne Rokita, Gregory P. Way
:contact: Jo Lynne Rokita (rokita@email.chop.edu)
:organization: CHOP
:status: complete
:date: |date|

.. meta::
   :keywords: pdx, mouse, WES, RNA-Seq, Fusions, SNP array, TMB, breakpoints, 2019
   :description: code to create PPTC PDX classifier plots and correlation figures

Introduction
============

Here, we provide scripts to enable reproducible generation of Manuscript Figures 4 and S4. This repo contains code for:

1. Creating ROC Curves from classifier data
2. Creating gene and variant level classification data
3. Creating correlation matrices for classification scores, tumor mutation burden, and mutational signatures.

Details
=======

- RUN-classifier-analyses.R
- classifier-plots-revision.R
- corplots.R
- cormat.R
- install.packages.R
- mutation-color-function.R
- theme.R




Software Requirements
=====================

R 3.4.3

Pipeline
========

.. code-block:: bash

         # How to run:
         # Download github repository in your home directory (~/)
         git clone https://github.com/marislab/pptc-pdx-classifier-analysis.git
   
         # Run script to create pie chart
         Rscript ~/create-pptc-pdx-oncoprints/R/RUN-classifier-analyses.R 
         
Results
========

Results will appear in the results folder

