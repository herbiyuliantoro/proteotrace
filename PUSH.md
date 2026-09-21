# Pushing this to GitHub

There is no `gh` CLI and no GitHub credentials on this machine, so the last
step has to run from your side. Two ways.

## Option A — with the GitHub CLI (easiest)

```sh
# install once
sudo apt install gh          # or: winget install GitHub.cli

gh auth login                # pick GitHub.com → HTTPS → browser

cd /mnt/c/ProteoTrace
gh repo create proteotrace --private --source=. --remote=origin --push
```

## Option B — without it

1. Create an **empty private** repo named `proteotrace` on github.com.
   Do not let it add a README, .gitignore or licence — this repo already has them.
2. Then:

```sh
cd /mnt/c/ProteoTrace
git remote add origin https://github.com/<your-username>/proteotrace.git
git push -u origin main
```

When it asks for a password, paste a personal access token, not your account
password. Create one at github.com → Settings → Developer settings → Personal
access tokens → Fine-grained, with **Contents: read and write** on this repo.

## After the first push

Two placeholders still say `<your-username>` / `<ORG-OR-USER>`:

```sh
cd /mnt/c/ProteoTrace
sed -i 's|<your-username>|YOURNAME|; s|<ORG-OR-USER>|YOURNAME|' README.md CITATION.cff
git commit -am "Point the repository links at the real GitHub account"
git push
```

## When the preprint goes out

1. Flip the repo to public.
2. Link it to Zenodo (zenodo.org → GitHub → toggle the repo on), then cut a
   release. Zenodo mints a DOI for that exact snapshot.
3. Put the bioRxiv DOI into `CITATION.cff` under `preferred-citation`.

That order matters. A public repo on its own gives you no timestamp anyone can
cite; the preprint DOI and the Zenodo DOI are what establish priority.
