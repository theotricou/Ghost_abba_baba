##### Appraoch 1 simulation #####

# git clone https://github.com/AADavin/Zombi.git

python3 ZOMBI/Zombi.py T Ghost_abba_baba/approach1/Parameters/zombi_parameters output_dir

python3 Ghost_abba_baba/approach1/build_ms_command.py output_dir/T/CompleteTree.nwk -p Ghost_abba_baba/approach1/Parameters/sim_parameters -s 20 -o output_dir -v

Rscript Ghost_abba_baba/approach1/ms_simulation.R output_dir
