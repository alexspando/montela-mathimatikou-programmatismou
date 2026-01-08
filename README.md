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

## Repository Structure (suggested)

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
