# Neutrino Induced Dimuon Event Generator

***

> Author : Sourav Sarkar

> Email : <ssarkar1@ualberta.ca>

> PhD Thesis Work : <https://github.com/ssarkarbht/PhDThesis>

***

## Generator Overview

***

This software framework is developed to run Monte Carlo (MC) based simulation of neutrino interactions producing two high energy outgoing muons from/near the interaction vertex, known as dimuon events. The two major physics processes producing these events are charm production in charge-current (CC) deep inelastic scattering (DIS) and trident interaction of muon neutrino netrinos ($\nu_{\mu}$). With the large statistics of detected neutrinos in neutrino telescopes, such as IceCube and P-ONE and in collider neutrino experiments such as FASER $\nu$ and SND@LHC, measurements on dimuon events can probe into various Standard Model and Beyond Standard Model tests, e.g. QCD Physics, Parton Density Functions, SM extension models, and SUSY models.

## Running the Event Generator

The following steps show how to run an example event generation chain for producing MC events. The full generation process can be broken down into multiple stages and the framework is highly modular and customizable for different needs of different experimental setups.

1. Get the singularity container

```
singularity pull --arch amd64 library://ssarkarbht/simgen/nudimuon-generator:v1.0
```

2. Launch the singularity shell with your directories mounted

```
singularity shell -B /data:/data -B <Other Paths>:<Other Paths> nudimuon-generator_v1.0.sif
```

This will take you to the container environment where you can run all the generator scripts.

3. Get the dimuon generator repo

```
git clone https://github.com/ssarkarbht/NuDimuon-Generator.git
```

4. Set up the generator paths and environments

```
cd NuDimuon-Generator
source setup.sh
```

5. Scripts for running individual step of the simulation process are located at `scripts/`. To run an end-to-end event generation process for either collider neutrino (e.g. FASERnu, SND@LHC) or neutrino telescope (e.g. IceCube, P-ONE) experiements, single bash scripts are available in `example/`. The bash script takes an input simulation configuration file (also located in the respective directories) which has the option to define all the tunable parameters of the event generation process.

6. For the mass production of the events in a remote cluster, example batch submission scripts for * HTCondor * workload manager is provided in the respective example subdirectories.

