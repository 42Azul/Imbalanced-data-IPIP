# Imbalanced-data-IPIP

This repository contains the implementation of an ensemble-based training algorithm for imbalanced data using the IPIP approach. The algorithm combines different training algorithms such as Support Vector Machines (SVM), Regularized Logistic Regression (RLOG), Gradient Boosting Machines (GBM), and Random Forest (RANGER) in R scripts. 

## Background

This project is a part of a Bachelor's thesis focused on handling imbalanced data in the context of health systems, specifically related to the COVID-19 dataset. The aim is to develop an effective method for training models with imbalanced data using the IPIP technique. While the current implementation utilizes SVM, RLOG, and RANGER, additional algorithms can be incorporated into the ensemble to further enhance performance.

## Methodology

The IPIP (Imbalanced Problem Identification and Prediction) approach is utilized to address the imbalanced data problem. The algorithm forms ensembles by combining multiple base learners, including SVM, RLOG, RANGER, and GBM. Each base learner brings its unique strengths and characteristics to the ensemble, allowing for a comprehensive analysis of the imbalanced data.

### Repository Structure

The repository is organized into the following folders:

- **Current Work**: This folder contains the latest developments and updates of the project. It includes the most recent code, documentation, and experiments.

- **Datasets**: The datasets folder contains the relevant data used for training and evaluating the models. It includes the COVID-19 dataset and any additional datasets used for comparison or analysis.

- **IPIP explanation**: This folder provides some previous explanations of the IPIP approach.

- **Related works**: In the Related works folder, you will find research papers and studies that explore the topic of imbalanced data. These resources provide valuable information, techniques, and findings from the broader scientific community, helping to shape the understanding of imbalanced data challenges and potential solutions.
