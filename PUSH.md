# Pushing this to GitHub

One thing to get straight first, because it causes most of the confusion:

**This folder has two names.** In Windows it is `C:\ProteoTrace`. In WSL the
same folder is `/mnt/c/ProteoTrace`. Same files, two shells. Pick one shell and
stay in it — Linux commands like `sudo apt` will not work in PowerShell, and
the `/mnt/c/...` path does not exist there.

The recommended route is WSL, because git is already set up there and nothing
needs installing.

---

## Route 1 — WSL, nothing to install ⭐

Git for Windows ships a credential manager, and this repo is already pointed at
it, so signing in happens in your browser. No token to create.

**Step 1.** Create the repository on GitHub, in your browser:
<https://github.com/new>

- Name: `proteotrace`
- Visibility: **Private**
- Leave "Add a README", ".gitignore" and "licence" **unticked** — this repo
  already has all three, and ticking them will cause a conflict.

**Step 2.** In a **WSL / Ubuntu** terminal:

```bash
cd /mnt/c/ProteoTrace
git remote add origin https://github.com/herbiyuliantoro/proteotrace.git
git push -u origin main
```

A browser window opens, you approve, and it pushes. Replace `herbiyuliantoro` with
your GitHub username.

---

## Route 2 — PowerShell with the GitHub CLI

If you would rather work in PowerShell, install the CLI first. Note `winget`,
not `apt`; `apt` is a Linux tool.

```powershell
winget install --id GitHub.cli
```

Close PowerShell and open a new one, so it picks up the new command. Then:

```powershell
cd C:\ProteoTrace
gh auth login
gh repo create proteotrace --private --source=. --remote=origin --push
```

This route creates the repository for you, so you can skip the browser step.

---

## After the first push

Two links still say `herbiyuliantoro`. In **WSL**:

```bash
cd /mnt/c/ProteoTrace
sed -i 's|herbiyuliantoro|herbiyuliantoro|; s|herbiyuliantoro|herbiyuliantoro|' README.md CITATION.cff
git commit -am "Point the repository links at the real GitHub account"
git push
```

---

## If the push is rejected

`Updates were rejected because the remote contains work that you do not have`
means GitHub added a README when the repository was created. Easiest fix is to
delete that repository and make a new empty one, following step 1 exactly.

---

## When the preprint goes out

1. Switch the repository to public.
2. Connect it to Zenodo: <https://zenodo.org/account/settings/github/>, toggle
   `proteotrace` on, then cut a release on GitHub. Zenodo mints a DOI for that
   exact snapshot.
3. Put the bioRxiv DOI into `CITATION.cff` under `preferred-citation`.

The order matters. A public repository by itself gives you nothing anyone can
cite. The preprint DOI and the Zenodo DOI are what establish priority.
