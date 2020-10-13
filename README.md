# wormjam-db
A repository containing metabolites and lipids from WormJam, in silico predictions and curated from literature

# Introduction

<i>Caenorhabditis elegans</i> (<i>C. elegans</i>) is an important model organism in biomedical research. Despite its importance no metabolome and/or lipidome exits. This repository is closing this gap by delivering metabolites curated from literature, predicted to be present or generated in silico.

# Structure

This repository contains different folder representing the different levels of curation.

## Literature

This folder contains metabolites curated from published <i>C. elegans</i> metabolomics and lipidomics studies. For each publication one file is provided with metabolites and names normalized towards a standard name (mostly ChEBI).

## Literature combined

The folder <code>literature_combined</code> contains a list of all unique metabolites across the entire curated literature.

## Prediction

This folder contains predicted metabolites, e.g. based on detected and identified metabolites, using variable scafolds such as acyl chains of different length.
