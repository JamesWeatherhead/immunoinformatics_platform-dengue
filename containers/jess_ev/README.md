# JessEV epitope selection container

This container runs the [JessEV](https://github.com/SchubertLab/JessEv) integer
linear program for epitope set selection, requiring a Gurobi optimizer license.

## Important: license credentials must be supplied at runtime

The upstream repo (`pmccaffrey6/immunoinformatics_platform`) inherited
committed Gurobi credentials in `src/gurobi.lic` and the `Dockerfile`'s `ENV`
variables. This fork (`JamesWeatherhead/immunoinformatics_platform-dengue`)
has scrubbed those credentials in HEAD; the upstream history still contains
the originals and those should be considered compromised and rotated by their
owner.

## How to supply credentials

### Option A: environment variables at `docker run` / `singularity exec`

```bash
docker run -e WLSACCESSID=$WLSACCESSID \
           -e WLSSECRET=$WLSSECRET \
           -e LICENSEID=$LICENSEID \
           jess_ev <command>
```

### Option B: mount your own license file

```bash
docker run -v /path/to/your/gurobi.lic:/gurobi.lic jess_ev <command>
```

### Option C: replace `src/gurobi.lic` with your own values before `docker build`

Edit `src/gurobi.lic` to substitute your real `WLSACCESSID`, `WLSSECRET`, and
`LICENSEID`. Do NOT commit the modified file (it is not in `.gitignore` here
because the placeholder values are themselves checked in; a developer note is
the only safeguard).

## Where to get a Gurobi license

Free academic licenses:
- https://www.gurobi.com/academia/academic-program-and-licenses/

Commercial WLS:
- https://www.gurobi.com/products/gurobi-optimizer/

## Original example invocations (preserved from upstream)

```bash
sudo docker run --rm jess_ev /bin/bash -c "/opt/conda/envs/jessev/bin/python3 /JessEV/design.py"

sudo docker run --rm -v /home/pathinformatics:/home/pathinformatics pmccaffrey6/jess_ev:latest \
  /bin/bash -c "/opt/conda/envs/jessev/bin/python3 /JessEV/design.py \
                /home/pathinformatics/jess_ev_csv.csv \
                /home/pathinformatics/output.csv"
```
