# Deploying AptaRank

The lab server (`aai`, Ubuntu 24.04, 256 cores) runs the dashboard so that
biologists can use AptaRank from a browser without installing anything. Being
Linux, it can also run fpocket — so Tier 2 works on real structures there,
which it cannot on Windows.

## What the deployment looks like

```
workstation                          server (aai)
───────────                          ────────────
git push ──────────────────────────► ~/aptarank.git      bare repo
                                     ~/aptarank          checkout at a known commit
                                       .venv/            its own Python environment
                                     ~/aptarank-data/    runs, uploads, caches, libraries
                                     ~/.local/bin/fpocket

browser ◄── ssh tunnel ────────────► 127.0.0.1:8501      Streamlit, loopback only
```

Three decisions worth knowing:

**Code comes from git; data never does.** The checkout is replaced on every
deploy, so it must contain nothing anyone would miss. Uploads, results and
reference libraries live in `~/aptarank-data`, which no deployment touches.

**The dashboard binds to loopback and is reached through an SSH tunnel.**
Streamlit has no authentication, and this service accepts uploads and launches
subprocesses; putting that on a network — even an internal one — with no login
is not defensible for unpublished sequence data. The tunnel makes the server's
existing SSH keys the authentication: per person, already managed, nothing new
to leak. It also encrypts the traffic. The cost is that each user needs a key
on the server; for a lab group that is the right trade.

**AptaRank takes a fixed share of a shared machine.** Two analyses at a time,
16 folding workers each — 32 of 256 cores — plus `nice 10` and no GPU access.
The limits live in `deploy/aptarank.sh` and override anything a user's settings
ask for.

## First deployment

From the workstation, with the server in `~/.ssh/config` (as `H200` here):

```powershell
.\deploy\deploy.ps1
```

That pushes the current commit, builds the environment (including ushuffle and
fpocket from source), runs the test suite, and starts the dashboard. It refuses
to start the service if the tests fail.

Then, to use it:

```powershell
.\deploy\connect.bat        # opens the tunnel and the browser
```

## Updating after a code change

```powershell
git commit -am "..."
.\deploy\deploy.ps1
```

The dashboard restarts. **Analyses already running are left alone** — they are
detached subprocesses that write their own results, and killing a half-finished
scientific run to restart a web page would be the wrong trade. They finish
under the code that started them, which is also why every artifact records its
own commit.

## Server commands

```bash
~/aptarank/deploy/aptarank.sh start|stop|restart|status|logs
```

`status` reports the deployed commit, the job count and whether fpocket is
available.

## Reference libraries and targets on the server

Put a validated corpus at `~/aptarank-data/data/corpus/<name>.csv`, with a
`<name>.manifest.json` beside it recording provenance:

```json
{
  "name": "Validated RNA aptamers v1",
  "source": "curated from SELEX literature",
  "curator": "Laura",
  "curated_date": "2026-08-01"
}
```

Without that manifest the library still works, but the dashboard marks it
"provenance not recorded" — correct columns are not provenance.

Build target evidence on the server, where fpocket is available:

```bash
cd ~/aptarank && .venv/bin/python -m aptarank target build \
    --pdb-id 3SPU --chain A \
    --set 'tier2.target.active_site_residues=[120,122,189]' \
    --set 'tier2.target.retain_hetero_resnames=["ZN"]' \
    -o ~/aptarank-data/cache/targets
.venv/bin/python scripts/verify_target_bundle.py ~/aptarank-data/cache/targets/3SPU_*.bundle.json
```

Review the output — chain, retained hetero groups, which cavity was selected
and why — before letting anyone use it. Bundle creation stays a deliberate step
rather than a button: a PDB ID alone does not settle the model, the chain, which
ligands to keep, or which cavity is the functional one, and one-click generation
would produce persuasive-looking evidence from decisions the biologist never
knew they made.

For numbers that go in the paper, prefer bundles from the pinned CI workflow
(`.github/workflows/target-bundle.yml`); the server is for exploration.

## If browser access without a tunnel becomes necessary

Do not simply bind Streamlit to `0.0.0.0`. Put an authenticating reverse proxy
in front of it:

1. Fetch the Caddy static binary into `~/.local/bin` (no root needed for ports
   above 1024).
2. Serve `https://<host>:8443` with per-person credentials — `caddy hash-password`
   for each, never one shared password.
3. Proxy to `127.0.0.1:8501` with WebSocket support (Caddy's `reverse_proxy`
   handles this by default) and a long read timeout, since analyses stream
   progress.
4. Use the institution's internal CA. Caddy's self-signed `tls internal` works
   only if its root certificate is installed as trusted on every client — an
   untrusted-certificate warning trains exactly the wrong habit.
5. Ask whoever administers the network to restrict the port to the lab subnet.

Until all five are in place, the tunnel is the safer option.

## Uninstalling

```bash
~/aptarank/deploy/aptarank.sh stop
rm -rf ~/aptarank ~/aptarank.git      # code
# ~/aptarank-data holds the results — delete it deliberately, if at all
```
