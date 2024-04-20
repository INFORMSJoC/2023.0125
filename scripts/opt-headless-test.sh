#!/bin/bash


# 01
jupyter nbconvert --to script opt-01-heur-sol-algo.ipynb
python3 opt-01-heur-sol-algo.py\
    --fmin=0 --fmax=5 --finc=1 --rhat=3\
    --casestudy="imelda" --approach="stochastic"\
    --numprocesses=8
python3 opt-01-heur-sol-algo.py\
    --fmin=0 --fmax=5 --finc=1 --rhat=3\
    --casestudy="imelda" --approach="robust"\
    --numprocesses=8
rm opt-01-heur-sol-algo.py


# 02
jupyter nbconvert --to script opt-02-heur-sol-eval.ipynb
for omega in $(seq 1 4); do
    python3 opt-02-heur-sol-eval.py\
        --fmin=0 --fmax=5 --finc=1 --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic" --omega=$omega\
        --numchunks=32 --numprocesses=8
done
for omega in $(seq 1 4); do
    python3 opt-02-heur-sol-eval.py\
        --fmin=0 --fmax=5 --finc=1 --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="robust" --omega=$omega\
        --numchunks=32 --numprocesses=8
done
rm opt-02-heur-sol-eval.py


# 03
jupyter nbconvert --to script opt-03-ws.ipynb
python3 opt-03-ws.py\
    --fmin=0 --fmax=5 --finc=1 --rhat=3\
    --casestudy="imelda" --pftype="dc"\
    --numprocesses=8
rm opt-03-ws.py


# 04
jupyter nbconvert --to script opt-04-ev.ipynb
for f in $(seq 0 5); do
    python3 opt-04-ev.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic"\
        --timelimit=600
done
rm opt-04-ev.py


# 05
jupyter nbconvert --to script opt-05-eev.ipynb
for f in $(seq 0 5); do
    python3 opt-05-eev.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic"\
        --timelimit=600
done
rm opt-05-eev.py


# 06
jupyter nbconvert --to script opt-06-mv.ipynb
for f in $(seq 0 5); do
    python3 opt-06-mv.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="robust"\
        --timelimit=600
done
rm opt-06-mv.py


# 07
jupyter nbconvert --to script opt-07-mmv.ipynb
for f in $(seq 0 5); do
    python3 opt-07-mmv.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="robust"\
        --timelimit=600
done
rm opt-07-mmv.py


# 08
jupyter nbconvert --to script opt-08-sp-ro.ipynb
for f in $(seq 0 5); do
    python3 opt-08-sp-ro.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic"\
        --timelimit=600
done
for f in $(seq 0 5); do
    python3 opt-08-sp-ro.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="robust"\
        --timelimit=600
done
rm opt-08-sp-ro.py


# 09
jupyter nbconvert --to script opt-09-no-good-cut.ipynb
for f in $(seq 1 5); do  # do not solve for f=0
    python3 opt-09-no-good-cut.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic"\
        --timelimit=600
done
for f in $(seq 1 5); do  # do not solve for f=0
    python3 opt-09-no-good-cut.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="robust"\
        --timelimit=600
done
rm opt-09-no-good-cut.py


# 10
jupyter nbconvert --to script opt-02-heur-sol-eval.ipynb
for omega in $(seq 1 4); do
    python3 opt-02-heur-sol-eval.py\
        --fmin=0 --fmax=5 --finc=1 --rhat=3\
        --pftype="lpacf" --casestudy="imelda" --approach="stochastic" --omega=$omega\
        --numchunks=32 --numprocesses=8
done
rm opt-02-heur-sol-eval.py
jupyter nbconvert --to script opt-08-sp-ro.ipynb
for f in $(seq 4 5); do
    python3 opt-08-sp-ro.py\
        --f=$f --rhat=3\
        --pftype="lpacf" --casestudy="imelda" --approach="stochastic"\
        --timelimit=900
done
rm opt-08-sp-ro.py
jupyter nbconvert --to script opt-10-cross-powerflow-eval.ipynb
python3 opt-10-cross-powerflow-eval.py\
    --fmin=0 --fmax=5 --finc=1 --rhat=3\
    --casestudy="imelda" --pftype-a="dc" --pftype-b="lpacf" --approach="stochastic"\
    --timelimit=900
rm opt-10-cross-powerflow-eval.py


# 11
jupyter nbconvert --to script opt-11-cross-uncertainty-eval.ipynb
python3 opt-11-cross-uncertainty-eval.py\
    --fmin=0 --fmax=5 --finc=1 --rhat=3\
    --casestudy="imelda" --pftype="dc"
rm opt-11-cross-uncertainty-eval.py


# 12
jupyter nbconvert --to script opt-12-acpf-eval.ipynb
for f in $(seq 0 5); do
    python3 opt-12-acpf-eval.py\
        --f=$f --rhat=3\
        --pftype="dc" --casestudy="imelda" --approach="stochastic"
done
rm opt-12-acpf-eval.py
