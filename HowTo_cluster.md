# HowTo cluster

Open a VPN-connection to the university network.

Connect to the cluster per ssh in the terminal:

```bash
ssh -X ostertms@ipvslogin.informatik.uni-stuttgart.de
ssh -X sgscl1
```

`ssh` connect with `Files` to my userspace:

```bash
ssh://ostertms@ipvslogin.informatik.uni-stuttgart.de
```

Run the job script with the filename of the scenario to execute:

```bash
./job.sh lid_driven_cavity
```

Copy the resulting `.vti` stack to the host computer.