# bumpy_GPU_CPP
Temporary stores the CUDA GPU implementation of the 2D bumpy compression code.

To compile in the interactive node use the following:

```console
srun --pty -p gpu -c 2 -t 2:00:00 --gres=gpu:v100:1 bash
module restore cuda10 
nvcc -x cu nested_for_loop.cpp -o test.exe -gencode arch=compute_70,code=sm_70
```
