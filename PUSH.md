# Releasing a new version

1. Set the same version in `pyproject.toml` and `proteotrace/__init__.py`.
2. Run the included example and regenerate the figures as described in the
   README. Check the validation instructions in `docs/TUTORIAL.md`.
3. Commit the changes, create a Git tag for that version, and publish a GitHub
   release from the tag.
4. Enable the repository in Zenodo if you want a DOI for the release.

Add a paper citation to `CITATION.cff` when a paper or preprint DOI is
available.
