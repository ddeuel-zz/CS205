# CS205
Grain Surface Chemistry

serial: gcc -DUSE_CLOCK grain.c timing.c -o grain -lm

multi-threaded CPU: pgcc -Minfo -ta=multicore -DUSE_CLOCK grain_acc_cpu.c timing.c -o grain_acc_cpu -lm

GPU: pgcc -acc -Minfo -DUSE_CLOCK grain_acc.c timing.c -o grain_acc -lm

To use OpenACC run these commands to configure the shell environment:

export PGI=/opt/pgi;

export PATH=/opt/pgi/linux86-64/19.4/bin:$PATH;

To change the number of threads used in CPU parallel execution set the environment variable like this:

export ACC_NUM_CORES=2 
