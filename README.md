[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Comparisons of Two-Stage Models for Flood Mitigation of Electrical Substations

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Comparisons of Two-Stage Models for Flood Mitigation of Electrical Substations](https://doi.org/10.1287/ijoc.2023.0125) by B. Austgen, E. Kutanoglu, J. J. Hasenbein, and S. Santoso. 
The snapshot is based on 
[this SHA](https://gitlab.com/BrentAustgen/two-stage-hurricane-model/-/commit/dabf5fa757810e69726b6bd96de985c5650b67a0)
in the development repository. 

**Important: This code is being developed on an on-going basis at 
https://gitlab.com/BrentAustgen/two-stage-hurricane-model. Please go there if you would like to
get a more recent version or would like support.**



## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0125

https://doi.org/10.1287/ijoc.2023.0125.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@misc{AustgenComparisons2024,
  author =        {B. Austgen and E. Kutanoglu and J. J. Hasenbein and S. Santoso},
  publisher =     {INFORMS Journal on Computing},
  title =         {Comparisons of Two-Stage Models for Flood Mitigation of Electrical Substations},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0125.cd},
  url =           {https://github.com/INFORMSJoC/2023.0125},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0125},
}  
```



## Description

The software in this repository allows users to the solve the optimization models found in
the paper and execute the various analyses we conducted on the results of from those optimization
models. The optimization models are designed to inform the siting and sizing of ad hoc flood
mitigation around electrical substations prior to an imminent hurricane such that the resulting
load shed may be minimized according to some preferred risk metric. The analyses mainly pertain
to the sensitivity of the models to various parameters and to the performance of the models'
prescribed mitigation solutions in an exact AC power flow model.



## Environment
The scripts we provide are compatible with the conda environment specified in `environment.yaml`,
though the file specifies numerous packages that are not called upon by this application.

Running the optimization scripts requires a working installation of Gurobi Optimizer and IPOPT.
The former is used to solve mixed-integer linear programs (MILPs) and
convex mixed-integer quadratically constrained programs (MIQPs) to global optimality.
The latter is used to solve a nonconvex optimal AC power flow model to local optimality.



## Data (`data/`)

To replicate some of the figures that appear in the paper, users of this repository must create a
directory called `data/original` and populate it with the following third-party data, which may
be downloaded as individual country files from 
[geoBoundaries @ William & Mary](https://www.geoboundaries.org/):
- `geoBoundaries-MEX-ADM1-all.zip` 
- `geoBoundaries-USA-ADM1-all.zip`

Files in `data/preprocessed` are those produced by our data preprocessing scripts, which are not
present in this repository. These files are used primarily to build instances of our two-stage
optimization models. They are also used to postprocess those results and to build some of the
figures that appear in the paper.



## Source (`src/`)

The source Python files contain helper classes and functions as well as filename definitions.
Below, we provide a brief description of each source file:
- `config.py`: directory and filename definitions for data and results
- `util.py`: tools for optimization and postprocessing
- `model.py`: implementations of the two-stage optimization models in the paper
- `drawutil.py`: tools to aid figure generation



## Scripts (`scripts/`)

We provide a variety of scripts in the form of Jupyter notebooks in the `scripts/` directory.
The purpose of each such script is detailed below.
- `opt-01-heur-sol-algo.ipynb`: run the greedy heuristic and stores the mitigation solutions
- `opt-02-heur-sol-eval.ipynb`: evaluate the greedy heuristic mitigation solutions in a power flow model
- `opt-03-ws.ipynb`: solve the expected wait-and-see (EWS) and maximum wait-and-see (MWS) models
- `opt-04-ev.ipynb`: compute an expected value (EV) solution
- `opt-05-eev.ipynb`: solve the expected result of the expected value solution (EEV) models
- `opt-06-mv.ipynb`: compute a maximum value (MV) solution
- `opt-07-mmv.ipynb`: solve the maximum result of the maximum value solution (MMV) models
- `opt-08-sp-ro.ipynb`: solve the stochastic programming (SP) and robust optimization (RO) models
- `opt-09-no-good-cut.ipynb`: solve the SP and RO models again after eliminating the optimal mitigation solution with a "no-good" cut to evaluate its uniqueness
- `opt-10-cross-powerflow-eval.ipynb`: evaluate performance of mitigation solution induced by a particular power flow model under an alternative model of power flow
- `opt-11-cross-uncertainty-eval.ipynb`: evaluate performance of SP mitigation solution in RO model and vice versa
- `opt-12-acpf-eval.ipynb`: evaluate the performance of the SP and RO solutions in the ACPF model
- `post-01-consolidate.ipynb`: consolidate certain results into CSV files
- `ijoc-01-lpac-geometries.ipynb`: illustrate the disc and cosine geometries of our LPAC variants
- `ijoc-02-objective-and-bounds.ipynb`: illustrate the EEV/SP/EWS MMV/RO/MWS objective value data
- `ijoc-03-solution-times.ipynb`: illustrate the time-to-solution for all SP and RO instances
- `ijoc-04-sp-ro-similarity.ipynb`: illustrate the similarity of the SP and RO mitigation solutions and the performance of one under the alternative uncertainty perspective
- `ijoc-05-uncertainty-heatmap.ipynb`: illustrate the flooding scenarios
- `ijoc-06-flood-and-mitigation-maps.ipynb`: illustrate a specific mitigation solution on a map
- `ijoc-07-acpf-perf.ipynb`: illustrate the performance of our models' prescribed mitigation solutions in the ACPF model and the relationship between contingency severity and error magnitude
- `ijoc-08-nomenclature.ipynb`: illustrate the impact of substation flooding on buses, transmission lines, and transformers


We provide an example shell script named `opt-headless-test.sh` that demonstrates how the
optimization scripts may be converted to Python scripts and executed from a command line interface.
Note that this shell script runs a small subset of all the experiments from our paper. As
detailed in the table below, we executed a large number of experiments, and it would not be
practical to execute all of these experiments in serial on a personal computing device. We
completed these experiments using SKX compute nodes from the Stampede2 cluster at the
Texas Advanced Computing Center (TACC), a system that prior to its retirement was fed computing
jobs from users via Simple Linux Utility for Resource Management (SLURM). To replicate all of our
results, we recommend using a comparable system and job scheduling tool.

We additionally provide shell scripts for running our postprocessing script and the
figure-generating scripts in headless mode. These shell scripts operate on data supplied in the
repository and are likely to terminate in seconds to minutes depending on your system. We name
these scripts `post-headless.sh` and `ijoc-headless.sh`.

*Note: The table below summarizes the two-stage models, power flow models, case studies, and budget
values that we evaluated in our experiments. For each pair of two-stage model and case study,
we evaluated for using each incorporated power flow model all integer budgets between zero and the
budget threshold (i.e., maximum useful budget) indicated in the two right-most columns. Budget
thresholds marked with an asterisk (\*) were determined via optimization whereas those without were
precomputed. In the row for EEV,* $\overline{\boldsymbol{x}}$ *denotes the EV solution we
identified. Likewise in the row for MMV,* $\widehat{\boldsymbol{x}}$ *denotes the MV solution we
identified.*


<table>
    <thead>
        <td></td>
        <td colspan="2" align="center"><b>Model</b></td>
        <td colspan="2" align="center"><b>Case Study</b></td>
    </thead>
    <thead>
        <td></td>
        <td align="center"><b>Name</b></td>
        <td align="center"><b>Formulation / Description</b></td>
        <td align="center"><b>Imelda</b></td>
        <td align="center"><b>Harvey</b></td>
    </thead>
    <tr>
        <td rowspan="6" align="right"><b>Two-Stage Model</b></td>
        <td align="center">EEV</td>
        <td>$$\sum_{\omega \in \Omega} \text{Pr}(\omega) \mathcal{L}(\overline{\boldsymbol{x}}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">9</td>
        <td align="right">66</td>
    </tr>
    <tr>
        <td align="center">SP</td>
        <td>$$\min_{\boldsymbol{x} \in \mathcal{X}} \sum_{\omega \in \Omega} \text{Pr}(\omega) \mathcal{L}(\boldsymbol{x}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">20</td>
        <td align="right">193</td>
    </tr>
    <tr>
        <td align="center">EWS</td>
        <td>$$\sum_{\omega \in \Omega} \text{Pr}(\omega) \min_{\boldsymbol{x} \in \mathcal{X}} \mathcal{L}(\boldsymbol{x}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">11</td>
        <td align="right">66</td>
    </tr>
    <tr>
        <td align="center">MMV</td>
        <td>$$\max_{\omega \in \Omega} \mathcal{L}(\widehat{\boldsymbol{x}}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">8</td>
        <td align="right">62</td>
    </tr>
    <tr>
        <td align="center">RO</td>
        <td>$$\min_{\boldsymbol{x} \in \mathcal{X}} \max_{\omega \in \Omega} \mathcal{L}(\boldsymbol{x}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">*9</td>
        <td align="right">*62</td>
    </tr>
    <tr>
        <td align="center">MWS</td>
        <td>$$\max_{\omega \in \Omega} \min_{\boldsymbol{x} \in \mathcal{X}} \mathcal{L}(\boldsymbol{x}, \boldsymbol{\xi}^\omega)$$</td>
        <td align="right">*5</td>
        <td align="right">*48</td>
    </tr>
    <tr>
        <td rowspan="4" align="right"><b>Power Flow Model</b></td>
        <td align="center">DC</td>
        <td align="center">sine is linear, cosine is 1, no reactive power<br>branch conductance is negligible</td>
        <td align="center">✓</td>
        <td align="center">✓</td>
    </tr>
    <tr>
        <td align="center">LPAC-C</td>
        <td align="center">sine is linear, cosine is 1,<br>4-sided polyhedral relaxation of unit circle</td>
        <td align="center">✓</td>
        <td align="center">✓</td>
    </tr>
    <tr>
        <td align="center">LPAC-F</td>
        <td align="center">sine is linear, 8-sided polyhedral relaxation of cosine,<br>12-sided polyhedral relaxation of unit circle</td>
        <td align="center">✓</td>
        <td align="center"></td>
    </tr>
    <tr>
        <td align="center">QPAC</td>
        <td align="center">sine is linear, quadratic relaxation of cosine,<br>exact quadratic representation of unit circle</td>
        <td align="center">✓</td>
        <td align="center"></td>
    </tr>
</table>



## Results (`results/`)

The subdirectories of `results/` contain result data as listed below.
- `heuristic/`: heuristic mitigation solutions and their performance in instances of the SP and RO models
- `ws/`: partial solutions to instances the wait-and-see (WS) model
- `ev/`: partial solutions to instances of the expected value (EV) model
- `eev/`: partial solutions to instances of the expected result of the expected value solution (EEV) model
- `mv/`: partial solutions to instances of the maximum value (MV) model
- `mmv/`: partial solutions to instances of the maximum result of the maximum value solution (MMV) model
- `sp/`: partial solutions to instances of the stochastic programming (SP) model
- `ro/`: partial solutions to instances of the robust optimization (RO) model
- `spsol-romod/`: evaluation of optimal mitigation solutions to instances of the SP model in corresponding RO instance
- `rosol-spmod/`: evaluation of optimal mitigation solutions to instances of the RO model in corresponding SP instance
- `acpf/`: evaluation of optimal mitigation solutions to instances of the SP and RO models under an exact AC power flow model
- `conslidated/`: condensed partial solution data in tabular format
- `figures/`: the JPG and EPS files that appear in the paper and its online supplement

Below, we describe and show all of the figures in `results/figures/` that appear in the paper and
its online supplement. We provide instructions for replicating each figure using the scripts and
data found in this repository.

Figure 1 in the paper shows the impact of substation flooding on buses and branches. It may be
replicated by running `scripts/ijoc-08-nomenclature.ipynb`.

![Figure 1](results/figures/ijoc-nomenclature.jpg)

Figure 2 in the paper shows the cosine and disc geometries for our three LPAC variants for the case
of $\beta_{nm} = 1.$ It may be replicated by running `scripts/ijoc-01-lpac-geometries.ipynb`.

![Figure 2](results/figures/ijoc-lpac-geometries.jpg)

Figure 3 in the paper shows the objective value and bounds as functions of the mitigation budget in
the Tropical Storm Imelda case study. It may be replicated by running
`scripts/ijoc-02-objective-and-bounds.ipynb`.

![Figure 3](results/figures/ijoc-obj-bnds-imelda.jpg)

Figure 4 in the paper shows the objective value and bounds as functions of the mitigation budget in
the Hurricane Harvey case study. It may be replicated by running
`scripts/ijoc-02-objective-and-bounds.ipynb`.

![Figure 4](results/figures/ijoc-obj-bnds-harvey.jpg)

Figure 5 in the paper shows the times required to solve the studied instances to optimality.
It may be replicated by running `scripts/ijoc-03-solution-times.ipynb`.

![Figure 5](results/figures/ijoc-solution-times.jpg)

Figure 6 in the paper shows the similarity of the SP and RO solutions and the performance of each
in the alternative model. It may be replicated by running `scripts/ijoc-04-sp-ro-similarity.ipynb`.

![Figure 6](results/figures/ijoc-sp-ro-similarity.jpg)

Figure 7 in the paper shows a comparison of the SP and RO mitigation solutions for the Hurricane
Harvey instance with $f=62$ and $\hat{r}=3$. It may be replicated by running
`scripts/ijoc-06-flood-and-mitigation-maps.ipynb` twice -- once specifying to illustrate the SP
solution and again specifying to illustrate the RO solution.

![Figure 7a](results/figures/ijoc-mitigation-harvey-dc-stochastic-r3-f62.jpg)
![Figure 7b](results/figures/ijoc-mitigation-harvey-dc-robust-r3-f62.jpg)

Figure 8 in the paper shows the absolute and relative error in load shed as a function of
contingency severity as determined by AC model evaluations of the optimal DC model solutions.
It may be replicated by running `scripts/ijoc-07-acpf-perf.ipynb`.

![Figure 8](results/figures/ijoc-acpf-error.jpg)

Figure 9 in the paper shows load shed resulting from the optimal DC model solutions evaluated
using DC and AC models. It may be replicated by running `scripts/ijoc-07-acpf-perf.ipynb`.

![Figure 9](results/figures/ijoc-acpf-perf.jpg)

Figure 1 in the online supplement shows a comparison of cosine relaxations based on tangent lines
that are equidistantly or optimally spaced. It may be replicated by running
`scripts/ijoc-01-lpac-geometries.ipynb`.

![Figure x](results/figures/ijoc-cosine-geometry-comparison.jpg)

Figure 2 in the online supplement shows the Tropical Storm Imelda flood levels by scenario and by
substation. It may be replicated by running `scripts/ijoc-05-uncertainty-heatmap.ipynb`.
![Figure x](results/figures/ijoc-imelda-uncertainty-heatmap.jpg)

Figure 3 in the online supplement shows the Hurricane Harvey flood levels by scenario and by
substation. It may be replicated by running `scripts/ijoc-05-uncertainty-heatmap.ipynb`.
![Figure x](results/figures/ijoc-harvey-uncertainty-heatmap.jpg)



## Ongoing Development

This code is being developed on an on-going basis at the author's
[GitLab site](https://gitlab.com/BrentAustgen/two-stage-hurricane-model).



## Support

For support in using this software, submit an
[issue](https://gitlab.com/BrentAustgen/two-stage-hurricane-model/-/issues/new).
