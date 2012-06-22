#!/bin/sh

K=0.75
temperature=0.05
min_dt=0.3
max_dt=0.75
number=10
time=1e11
seed=2342

echo "Running simulation (random seed = $seed)..."
time ./simul --K $K --temperature $temperature \
	     --min-dt $min_dt --max-dt $max_dt --number $number \
	     --time $time --seed $seed
echo "Simulation finished."
xmessage -center "Simulation finished."
