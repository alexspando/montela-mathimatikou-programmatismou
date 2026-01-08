# Power System Expansion Planning — Decomposition Methods (Julia/JuMP)

Semester project repository for the course **Mathematical Models Programming (9th semester)**,  
School of Electrical & Computer Engineering, **National Technical University of Athens (NTUA)**.

This project implements and documents classical decomposition algorithms for **Power System Expansion Planning (PSEP)** in deterministic and stochastic settings.

---

## Overview

The goal is to solve the expansion planning problem using:

- **Benders Decomposition** (deterministic case)
- **L-shaped Method** (two-stage stochastic programming)
- **Multi-cut L-shaped Method** (scenario-wise cuts for faster convergence)

All implementations are written in **Julia** using **JuMP** and **GLPK**, and produce `.csv` outputs for the iteration history and convergence tracking.

---

## Repository Structure 

```text
.
├── report.tex
├── report.pdf
├── code/
│   ├── benders.jl
│   ├── lshaped.jl
│   └── lshaped_multicut.jl
├── data/
│   ├── technology.csv
│   └── needs.csv
├── results/
│   ├── benders_results.csv
│   ├── lshaped_results.csv
│   └── lshaped_multicut_results.csv
└── README.md
Requirements

Julia ≥ 1.8

Packages:

JuMP

GLPK

CSV

DataFrames

Printf

Install dependencies:

using Pkg
Pkg.add(["JuMP", "GLPK", "CSV", "DataFrames"])

Data Files

The code expects the following input files:

data/technology.csv
Contains technology names, marginal costs, and investment costs.

data/needs.csv
Contains the load-duration curve slicing (durations and load levels).

Make sure the CSV header names match those referenced in the Julia scripts.

How to Run

Clone the repository:

git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>


Run each method from Julia (examples):

Benders Decomposition
include("code/benders.jl")


Outputs:

results/benders_results.csv (iteration log and convergence)

L-shaped Method
include("code/lshaped.jl")


Outputs:

results/lshaped_results.csv

Multi-cut L-shaped Method
include("code/lshaped_multicut.jl")
