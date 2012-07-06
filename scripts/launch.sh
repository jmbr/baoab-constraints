#!/bin/sh

temperature=0.5
friction=1e4
min_dt=0.08
max_dt=0.15
number=10
equilibration_time=1e5
production_time=1e11
seed=$RANDOM

echo "$seed" > random-seed.dat  # Save random seed.

# Run simulation.
echo "Running simulation (random seed = $seed)..."
time ./simul --temperature $temperature --friction $friction \
             --min-dt $min_dt --max-dt $max_dt --number $number \
             --equilibration $equilibration_time --time $production_time \
             --seed $seed

# Notify the user.
echo "Simulation finished."
xmessage -center "Simulation finished."
