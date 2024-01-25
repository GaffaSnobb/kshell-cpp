# Install
## ROCm
Tested with:
- `Ubuntu 22.04.3` kernel `6.5.0-15-generic` (however, AMD's official support with `22.04.3` is `6.2`).
- AMD RX 7900XTX, AMD 7950X, B650 chipset.

Follow the official ROCm installation instruction at: https://rocm.docs.amd.com/projects/install-on-linux/en/latest/index.html. Use the `amdgpu-install` route (instead of using `apt`). I installed the following usecases:
```
amdgpu-install --usecase=rocm,rocmdevtools,lrt,opencl,hip,hiplibsdk,openmpsdk
```
All of them might not be needed, I'll update the requirements when I have tested more. I had to install
```
apt install libstdc++-12-dev
```
for `hipcc` to find `cmath`. The exact version of `libstdc++` might be different for other situations, more info: https://github.com/ROCm/ROCm/issues/1843.
